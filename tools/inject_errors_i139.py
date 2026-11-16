"""I139 prereg 3.2 step 5: controlled error injection (60 errors) and detection.
Usage: python -I tools/inject_errors_i139.py
- 15 papers from extraction/sample_47.tsv, ordered by sha256("20261114-i139-inject:" + rec_id).
- Per paper 2 crude + 2 plausible errors; each error goes into its own copy of the reconciled JSON
  (extraction/inject/<rec_id>__<k>.json), so detection is attributed to one error.
- Crude: wrong estimate / CI bound / N / cycle / weight variable; the real quote is kept (it does not contain the planted value).
- Plausible: (P1) estimate and model label of another model in the same table, with a real text line that contains it;
  (P2) covariate list of another model; (P3) drop one confounder-concept covariate that the methods section mentions.
- Detection: (a) the step-3 checker reports a failure the clean copy does not have; or (b) the planted field disagrees
  with the independent extraction pass (pass B, primary; pass A, sensitivity). Covariates compared at concept level
  (tools/concepts_i139.py, prereg 5.4), numbers within 0.005, N within 0.5, cycles/weight via the compare normaliser.
Writes extraction/inject/injection_log.tsv and prints rates with Wilson 95% CIs.
"""
import sys, os, json, csv, hashlib, copy, re, math, random, subprocess
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from concepts_i139 import concept_set
from compare_extraction_i139 import norm as cnorm, get

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = "20261114-i139-inject"
OUT = os.path.join(ROOT, "extraction", "inject")
G_CONCEPTS = ["race", "education", "pir", "smoking", "alcohol", "physical_activity", "bmi_waist", "diabetes_dx", "lipid_drugs", "hypertension_dx"]


def wilson(k, n):
    if n == 0:
        return (0, 0)
    z = 1.959963984540054
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def t4(x):
    return int(x * 10000) / 10000


def fmt(v):
    return ("%.2f" % v) if abs(v) < 100 else ("%.1f" % v)


def text_lines(rid):
    return open(os.path.join(ROOT, "fulltext_txt", rid + ".txt"), encoding="utf-8", errors="replace").read().split("\n")


def find_line(lines, *needles):
    for ln in lines:
        if all(n in ln for n in needles) and len(ln) < 600:
            return ln.strip()
    return None


def checker_fails(path, rid):
    r = subprocess.run([sys.executable, "-I", os.path.join(ROOT, "tools", "check_extraction_i139.py"),
                        os.path.join(ROOT, "fulltext_txt", rid + ".txt"), path], capture_output=True, text=True, encoding="utf-8")
    return {ln.strip() for ln in r.stdout.splitlines() if ln.startswith("   ")}


def crude_options(j, rng):
    me = j["main_estimate"]
    opts = []
    est = me["estimate"]["value"]
    if isinstance(est, (int, float)):
        new = round(est * 1.18 + 0.03, 2) if est > 0 else round(est - 0.2, 2)
        opts.append(("C_estimate", [("main_estimate.estimate", new)]))
        lo = me["ci_low"]["value"]
        if isinstance(lo, (int, float)):
            opts.append(("C_ci_low", [("main_estimate.ci_low", round(lo * 0.85 - 0.02, 2))]))
    n = j["n_analytic"]["value"]
    if isinstance(n, int):
        opts.append(("C_n", [("n_analytic", int(n * 1.12) + 7)]))
    cyc = j["cycles"]["value"]
    if cyc:
        yrs = [int(y) for y in re.findall(r"(?:19|20)\d\d", " ".join(map(str, cyc)))]
        lo_y, hi_y = min(yrs), max(yrs)
        new = [f"{y}-{y + 1}" for y in range(lo_y + 2, hi_y + 2, 2)]  # drop the first cycle, add one after the last
        opts.append(("C_cycles", [("cycles", new)]))
    w = str(j["weight_variable"]["value"] or "")
    new_w = "WTINT2YR" if "MEC" in w.upper() else "WTSAF2YR"
    opts.append(("C_weight", [("weight_variable", new_w)]))
    rng.shuffle(opts)
    return opts


def plausible_options(j, rid, rng):
    me = j["main_estimate"]
    lines = text_lines(rid)
    opts = []
    est = me["estimate"]["value"]
    for o in j.get("other_models") or []:
        e = o.get("estimate")
        if isinstance(e, (int, float)) and isinstance(est, (int, float)) and abs(e - est) > 0.005 and o.get("ci_low") is not None:
            s_e = fmt(e)
            q = me["estimate"]["quote"] if s_e in (me["estimate"]["quote"] or "") else find_line(lines, s_e, fmt(o["ci_low"]))
            if q:
                lab_q = find_line(lines, str(o.get("model_label"))) or q
                opts.append(("P_other_model", [("main_estimate.estimate", e, q), ("main_estimate.ci_low", o["ci_low"], q),
                                               ("main_estimate.ci_high", o.get("ci_high"), q),
                                               ("main_estimate.model_label", o.get("model_label"), lab_q)]))
                break
    base = concept_set(j.get("covariates_final"))
    for o in j.get("other_models") or []:
        cv = o.get("covariates") or []
        if cv and concept_set(cv) != base:
            opts.append(("P_other_cov", [("covariates_final", list(cv))]))
            break
    meth = concept_set(get(j, "covariates_methods") or [])
    fin = j.get("covariates_final") or []
    droppable = [c for c in G_CONCEPTS if c in base and c in meth]
    rng.shuffle(droppable)
    for c in droppable[:2]:
        kept = [x for x in fin if c not in concept_set([x])]
        if concept_set(kept) != base:
            opts.append(("P_drop_cov", [("covariates_final", kept)], c))
    return opts


def set_field(j, path, value, quote=None):
    if path == "covariates_final":
        j["covariates_final"] = value
        return
    node = j
    for k in path.split(".")[:-1]:
        node = node[k]
    f = node[path.split(".")[-1]]
    f["value"] = value
    if quote is not None:
        f["quote"] = quote


def same(path, a, b):
    if path == "covariates_final":
        return concept_set(a) == concept_set(b)
    if path in ("main_estimate.estimate", "main_estimate.ci_low", "main_estimate.ci_high"):
        try:
            return abs(float(a) - float(b)) <= 0.005
        except Exception:
            return False
    if path == "n_analytic":
        try:
            return abs(int(a) - int(b)) < 1
        except Exception:
            return False
    return cnorm(a, path) == cnorm(b, path)


def val(j, path):
    return j.get("covariates_final") if path == "covariates_final" else get(j, path)


def main():
    os.makedirs(OUT, exist_ok=True)
    S = [r["rec_id"] for r in csv.DictReader(open(os.path.join(ROOT, "extraction", "sample_47.tsv"), encoding="utf-8"), delimiter="\t")]
    order = sorted(S, key=lambda r: hashlib.sha256(f"{SEED}:{r}".encode()).hexdigest())
    rows, papers = [], []
    for rid in order:
        if len(papers) == 15:
            break
        J = json.load(open(os.path.join(ROOT, "extraction", "pass_R", rid + ".json"), encoding="utf-8"))
        rng = random.Random(f"{SEED}:{rid}")
        cr = crude_options(J, rng)[:2]
        pl = plausible_options(J, rid, rng)
        kinds, chosen = set(), []
        for o in pl:  # prefer two different plausible kinds
            if o[0] not in kinds:
                chosen.append(o)
                kinds.add(o[0])
        for o in pl:
            if o not in chosen:
                chosen.append(o)
        pl = chosen[:2]
        if len(cr) < 2 or len(pl) < 2:
            print(f"skip {rid}: crude {len(cr)} plausible {len(pl)} (logged as deviation)")
            continue
        papers.append(rid)
        A = json.load(open(os.path.join(ROOT, "extraction", "pass_A", rid + ".json"), encoding="utf-8"))
        B = json.load(open(os.path.join(ROOT, "extraction", "pass_B", rid + ".json"), encoding="utf-8"))
        clean_path = os.path.join(ROOT, "extraction", "pass_R", rid + ".json")
        clean_fails = checker_fails(clean_path, rid)
        for k, (grp, opt) in enumerate([("crude", o) for o in cr] + [("plausible", o) for o in pl]):
            kind, edits = opt[0], opt[1]
            j = copy.deepcopy(J)
            for e in edits:
                set_field(j, e[0], e[1], e[2] if len(e) > 2 else None)
            p = os.path.join(OUT, f"{rid}__{k}.json")
            json.dump(j, open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
            new_fails = checker_fails(p, rid) - clean_fails
            fields = [e[0] for e in edits]
            disB = any(not same(f, val(j, f), val(B, f)) for f in fields)
            disA = any(not same(f, val(j, f), val(A, f)) for f in fields)
            base_agreeB = all(same(f, val(J, f), val(B, f)) for f in fields)
            rows.append({"rec_id": rid, "copy": k, "group": grp, "kind": kind, "fields": ",".join(fields),
                         "planted": json.dumps([e[1] for e in edits], ensure_ascii=False)[:200],
                         "dropped_concept": opt[2] if len(opt) > 2 else "",
                         "checker_new": " ; ".join(sorted(new_fails))[:300], "checker": int(bool(new_fails)),
                         "disagree_B": int(disB), "disagree_A": int(disA), "baseline_R_eq_B": int(base_agreeB),
                         "detected_B": int(bool(new_fails) or disB), "detected_A": int(bool(new_fails) or disA)})
    fl = open(os.path.join(OUT, "injection_log.tsv"), "w", encoding="utf-8", newline="")
    w = csv.DictWriter(fl, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    fl.close()
    print("papers:", " ".join(papers))
    for lab, key in (("primary (checker or pass B)", "detected_B"), ("sensitivity (checker or pass A)", "detected_A"), ("checker only", "checker")):
        for g in ("all", "crude", "plausible"):
            sub = [r for r in rows if g == "all" or r["group"] == g]
            k = sum(r[key] for r in sub)
            lo, hi = wilson(k, len(sub))
            print(f"{lab:34s} {g:9s} {k}/{len(sub)} = {t4(k / len(sub))}  Wilson95 [{t4(lo)}, {t4(hi)}]")
    sub = [r for r in rows if r["baseline_R_eq_B"]]
    k = sum(r["detected_B"] for r in sub)
    lo, hi = wilson(k, len(sub))
    print(f"conditional on clean R = pass B on the planted fields: {k}/{len(sub)} = {t4(k / len(sub)) if sub else 0}  Wilson95 [{t4(lo)}, {t4(hi)}]")
    for r in rows:
        if not r["detected_B"]:
            print("MISSED", r["rec_id"], r["group"], r["kind"], r["fields"], r["planted"], r["dropped_concept"])


if __name__ == "__main__":
    main()

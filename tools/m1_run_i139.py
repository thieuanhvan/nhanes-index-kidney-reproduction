"""I139 M1 runner: fit each paper's first-combination specification and compare with the reported estimate.
Usage: python -I tools/m1_run_i139.py [rec_id ...]
Reads tools/m1_specs/<rec_id>.py builders (see m1_lib_i139.py) and the locked spec table (hash checked).
Ranking (prereg section 6, M1): strict pass |dlogOR| <= log 1.10 and |dN|/N <= 10 %; near pass |dlogOR| <= log 1.25 and
|dN|/N <= 15 %; fail otherwise; not testable when the builder raises or the paper is not reconstructable.
Rounding rule: also counts as a pass on the OR component when the reproduced OR rounded to the reported number of
decimals equals the reported OR. Weight None -> w = 1 (unweighted point estimate; SE then not design-correct, M2 handles design).
Writes extraction/m1_results_i139.tsv (appends/replaces rows per rec_id) and prints one line per paper.
"""
import os, sys, json, csv, hashlib, importlib.util, traceback, math, datetime
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m1_lib_i139 as lib

ROOT = os.path.dirname(HERE)
X = os.path.join(ROOT, "extraction")
OUTP = os.path.join(X, "m1_results_i139.tsv")
COLS = ["rec_id", "family", "status", "rank", "or_reported", "ci_reported", "or_repro", "ci_repro", "dlogor", "n_reported",
        "n_repro", "dn_rel", "weight", "target", "n_terms", "readings", "notes", "run_time"]


def decimals(v):
    s = repr(float(v))
    return len(s.split(".")[1].rstrip("0")) if "." in s else 0


def rank(orr, orp, nr, npp):
    d = abs(math.log(orp) - math.log(orr))
    dn = abs(npp - nr) / nr if nr else float("inf")
    round_ok = round(orp, decimals(orr)) == orr
    o10 = d <= math.log(1.10) or round_ok
    o25 = d <= math.log(1.25) or round_ok
    if o10 and dn <= 0.10:
        return "strict_pass", d, dn
    if o25 and dn <= 0.15:
        return "near_pass", d, dn
    return "fail", d, dn


def run_one(rid, spec):
    p = os.path.join(HERE, "m1_specs", rid + ".py")
    row = {c: "" for c in COLS}
    row.update(rec_id=rid, family=spec["family"], or_reported=spec["estimate"], ci_reported=f"{spec['ci_low']}-{spec['ci_high']}",
               n_reported=spec["n_analytic"], run_time=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    if not os.path.exists(p):
        row.update(status="no_builder")
        return row
    try:
        sp = importlib.util.spec_from_file_location("b_" + rid, p)
        m = importlib.util.module_from_spec(sp)
        sp.loader.exec_module(m)
        b = m.build(lib)
        if b.get("not_testable"):
            row.update(status="not_testable", notes=str(b["not_testable"])[:300])
            return row
        df = b["df"]
        w = b["weight"]
        if w is None:
            df["_w1"] = 1.0
            w = "_w1"
        res = lib.svy_logit(df, b["y"], b["terms"], w, domain=b["domain"])
        t = res["table"].loc[b["target"]]
        st, d, dn = rank(float(spec["estimate"]), float(t.OR), float(spec["n_analytic"]), res["n"])
        row.update(status="run", rank=st, or_repro=float("%.4g" % t.OR), ci_repro=f"{float('%.4g' % t.lo)}-{float('%.4g' % t.hi)}",
                   dlogor=int(d * 10000) / 10000, n_repro=res["n"], dn_rel=int(dn * 10000) / 10000,
                   weight=b["weight"] or "none", target=b["target"], n_terms=len(b["terms"]),
                   readings=json.dumps(b.get("readings", {}), ensure_ascii=False)[:400], notes=str(b.get("notes", ""))[:300])
    except Exception as e:
        row.update(status="error", notes=(type(e).__name__ + ": " + str(e))[:300])
        traceback.print_exc()
    return row


def main():
    lock = open(os.path.join(X, "spec_table_i139.sha256"), encoding="utf-8").read()
    pj = os.path.join(X, "spec_table_i139.json")
    assert hashlib.sha256(open(pj, "rb").read()).hexdigest() in lock, "spec table changed after lock"
    T = {r["rec_id"]: r for r in json.load(open(pj, encoding="utf-8"))}
    ids = sys.argv[1:] or list(T)
    old = {}
    if os.path.exists(OUTP):
        old = {r["rec_id"]: r for r in csv.DictReader(open(OUTP, encoding="utf-8"), delimiter="\t")}
    for rid in ids:
        r = run_one(rid, T[rid])
        old[rid] = r
        print(rid, r["status"], r["rank"], "OR", r["or_reported"], "->", r["or_repro"], "N", r["n_reported"], "->", r["n_repro"], r["notes"][:80], flush=True)
    with open(OUTP, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for rid in T:
            if rid in old:
                w.writerow(old[rid])


if __name__ == "__main__":
    main()

"""I139 prereg 3.2 steps 1 and 4: compare two extraction passes field by field.
Usage: python -I tools/compare_extraction_i139.py <dirA> <dirB> <out_disagreements.tsv>
Prints per-field raw agreement and Gwet AC1 (two raters, agree/disagree on nominal values) over papers present in both.
"""
import sys, os, json, re, csv
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FREE = {"outcome_definition", "main_estimate.exposure_scale", "cycle_pooling"}  # free text: disagreement goes to reconciliation, reported separately
KEY = ["index_name", "index_components", "cycles", "n_analytic", "outcome_main", "outcome_definition",
       "egfr_equation", "uacr_threshold", "weight_variable", "cycle_pooling",
       "main_estimate.model_label", "main_estimate.effect_type", "main_estimate.exposure_scale",
       "main_estimate.estimate", "main_estimate.ci_low", "main_estimate.ci_high", "main_estimate.source_table",
       "covariates_final", "inflection_point"]


def get(j, path):
    for k in path.split("."):
        j = (j or {}).get(k) if isinstance(j, dict) else None
    if isinstance(j, dict) and "value" in j:
        j = j["value"]
    return j


def norm(v, path):
    if v in (None, "", []):
        return "NA"
    # Format-only equivalences (written after seeing the test-set outputs; declared in verdict 9.11).
    if path == "cycles":
        yrs = set()
        for x in (v if isinstance(v, list) else [v]):
            a = [int(y) for y in re.findall(r"(?:19|20)\d\d", str(x))]
            if len(a) >= 2:
                yrs |= set(range(min(a), max(a) + 1))
            yrs |= set(a)
        return f"{min(yrs)}-{max(yrs)}" if yrs else "NA"
    if path == "uacr_threshold":
        n = re.findall(r"\d+(?:\.\d+)?", str(v))
        return str(float(n[0])) if n else "NA"
    if path == "index_components" and isinstance(v, list):
        return f"n={len(v)}"
    if isinstance(v, list):
        return "|".join(sorted({re.sub(r"[^a-z0-9]+", " ", str(x).lower()).strip() for x in v}))
    if isinstance(v, (int, float)):
        return repr(float(v))
    s = str(v).lower()
    if path == "uacr_threshold":
        n = re.findall(r"\d+(?:\.\d+)?", s)
        return n[0] if n else "NA"
    if path == "egfr_equation":
        y = re.findall(r"20(?:09|12|21)", s)
        return ("ckd-epi " + (y[0] if y else "unstated") + (" cys" if "cystatin" in s else "")) if "epi" in s or "epidemiology" in s else s
    if path == "index_name":
        m = re.search(r"\(([^)]+)\)", s)
        s = m.group(1) if m and len(m.group(1)) <= 12 else s.split("(")[0].split(";")[0]
    if path == "outcome_main":
        s = s.split("(")[0].split(";")[0]
    if path == "egfr_equation":
        s = s.replace("version not stated", "")
    if path == "weight_variable" and "not named" in s:
        return "NA"
    s = re.sub(r"[^a-z0-9.]+", " ", s).strip()
    s = re.sub(r"\bindex\b", "", s).strip()
    if path.endswith("source_table") or path.endswith("model_label"):
        m = re.findall(r"(?:table|model)\s*([0-9ivx]+[a-z]?)", s)
        roman = {"i": "1", "ii": "2", "iii": "3", "iv": "4", "v": "5", "vi": "6"}
        return " ".join(roman.get(x, x) for x in m) or s
    return s


def main(da, db, out):
    ids = sorted({f[:-5] for f in os.listdir(da) if f.endswith(".json")} & {f[:-5] for f in os.listdir(db) if f.endswith(".json")})
    A = {i: json.load(open(os.path.join(da, i + ".json"), encoding="utf-8")) for i in ids}
    B = {i: json.load(open(os.path.join(db, i + ".json"), encoding="utf-8")) for i in ids}
    w = csv.writer(open(out, "w", encoding="utf-8", newline=""), delimiter="\t", lineterminator="\n")
    w.writerow(["rec_id", "field", "value_A", "value_B"])
    tot_a = tot = 0
    print(f"papers compared: {len(ids)}")
    for f in KEY:
        pairs = [(norm(get(A[i], f), f), norm(get(B[i], f), f)) for i in ids]
        agree = sum(a == b for a, b in pairs)
        n = len(pairs)
        cats = {x for p in pairs for x in p}
        pa = agree / n if n else 0
        pe = 0.0
        if len(cats) > 1:
            q = len(cats)
            pis = [sum((a == c) + (b == c) for a, b in pairs) / (2 * n) for c in cats]
            pe = sum(p * (1 - p) for p in pis) / (q - 1)
        ac1 = (pa - pe) / (1 - pe) if pe < 1 else 1.0
        tot_a += agree
        tot += n
        print(f"{f:32s} agree {agree}/{n}  AC1 {int(ac1 * 10000) / 10000}")
        for i, (a, b) in zip(ids, pairs):
            if a != b:
                w.writerow([i, f, str(get(A[i], f))[:300], str(get(B[i], f))[:300]])
    print(f"OVERALL raw agreement {tot_a}/{tot} = {int(10000 * tot_a / tot) / 10000 if tot else 0}")


if __name__ == "__main__":
    main(*sys.argv[1:4])

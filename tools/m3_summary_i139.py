"""I139 M3 summary: E3a(A) primary, E3a(B), secondary denominators, controls (rules 8-11 of extraction/m2_m3_rules_i139.md).
Usage: python -I tools/m3_summary_i139.py   -> prints and writes extraction/m3_summary_i139.json
Stratum weights N_h* = N_h x full-text eligibility proportion (verdict 9.9/9.11): N_h F1 7, F2 30, F3 49, F4 43, F5 17;
full-text eligible / assessed: F1 4/4, F2 13/13, F3 11/12, F4 11/13, F5 12/13.
"""
import os, sys, json, glob, math, hashlib
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C

NH = {"F1_tyg": 7, "F2_inflam": 30, "F3_diet_life": 49, "F4_lipid_anthro": 43, "F5_other": 17}
ELIG = {"F1_tyg": (4, 4), "F2_inflam": (13, 13), "F3_diet_life": (11, 12), "F4_lipid_anthro": (11, 13), "F5_other": (12, 13)}
NHS = {h: NH[h] * ELIG[h][0] / ELIG[h][1] for h in NH}
Z = 1.959963984540054
tr = lambda v: math.floor(v * 1e4) / 1e4 if v is not None and np.isfinite(v) else v


def wilson(k, n):
    if n == 0:
        return None
    p = k / n; d = 1 + Z * Z / n; c = (p + Z * Z / (2 * n)) / d; h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return [k, n, tr(p), tr(max(0.0, c - h)), tr(min(1.0, c + h))]


def stratified(items, seed="20261114-i139-e3a-boot"):
    """items: list of (stratum, 0/1). Weighted proportion with fpc, logit CI; bootstrap if a stratum is 0 or 1."""
    by = {}
    for h, v in items:
        by.setdefault(h, []).append(v)
    tot = sum(NHS[h] for h in by)
    W = {h: NHS[h] / tot for h in by}
    p = sum(W[h] * np.mean(by[h]) for h in by)
    var = sum(W[h] ** 2 * (1 - len(by[h]) / NHS[h]) * np.var(by[h], ddof=1) / len(by[h]) for h in by if len(by[h]) > 1)
    extreme = any(np.mean(by[h]) in (0.0, 1.0) for h in by)
    if not extreme and 0 < p < 1:
        se = math.sqrt(var); lg = math.log(p / (1 - p)); s = se / (p * (1 - p))
        lo, hi, method = 1 / (1 + math.exp(-(lg - Z * s))), 1 / (1 + math.exp(-(lg + Z * s))), "logit"
    else:
        rng = np.random.default_rng(int(hashlib.sha256(seed.encode()).hexdigest()[:15], 16))
        bs = []
        for _ in range(2000):
            bs.append(sum(W[h] * np.mean(rng.choice(by[h], len(by[h]), replace=True)) for h in by))
        lo, hi, method = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), "stratified bootstrap 2000"
    return {"weighted": tr(p), "lo": tr(lo), "hi": tr(hi), "method": method,
            "by_stratum": {h: [int(sum(by[h])), len(by[h])] for h in sorted(by)}}


def main():
    S = C.spec_table()
    R = {}
    for f in glob.glob(os.path.join(C.X, "m3", "*.json")):
        d = json.load(open(f, encoding="utf-8"))
        R[d["rec_id"]] = d
    L = C.labels()
    run = [d for d in R.values() if d["status"] == "run"]
    err = [d["rec_id"] for d in R.values() if d["status"] == "error"]
    out = {"n_files": len(R), "n_run": len(run), "errors": err}
    fp = lambda d, k: d[k]["class"] == "firm_positive"
    pp = lambda d, k: d[k]["class"] in ("firm_positive", "point_positive")
    unfl = [d for d in run if not d["flag_outcome_component"]]
    fl = [d for d in run if d["flag_outcome_component"]]
    for dag in ("A", "B"):
        k = "E3a_" + dag
        out[k + "_primary"] = stratified([(d["family"], int(fp(d, k))) for d in unfl])
        out[k + "_primary_unweighted"] = wilson(sum(fp(d, k) for d in unfl), len(unfl))
        out[k + "_point_estimate"] = stratified([(d["family"], int(pp(d, k))) for d in unfl])
        out[k + "_flagged"] = wilson(sum(fp(d, k) for d in fl), len(fl))
        near = [d for d in run if d["rec_id"] in NEAR]
        out[k + "_near_pass"] = wilson(sum(fp(d, k) for d in near), len(near))
        # whole sample: not testable with a missing label counted negative (lower) / positive (upper)
        nt = [r for r in S if r not in R or R[r]["status"] != "run"]
        miss_nt = [r for r in nt if L[r]["label_" + dag] == "missing"]
        k_run = sum(fp(d, k) for d in run)
        out[k + "_whole_sample_bounds"] = {"lower": wilson(k_run, len(S)), "upper": wilson(k_run + len(miss_nt), len(S)),
                                            "not_testable_missing": miss_nt}
        out[k + "_classes"] = {c: sorted(d["rec_id"] for d in run if d[k]["class"] == c)
                               for c in ("firm_positive", "point_positive", "negative")}
        out[k + "_deltas_run"] = {d["rec_id"]: [tr(d[k]["delta"]), tr(d[k]["lo"]), tr(d[k]["hi"])] for d in run if "delta" in d[k]}
    out["label_changes_A_to_B"] = sorted(d["rec_id"] for d in run if d["E3a_A"]["class"] != d["E3a_B"]["class"])
    out["neg_control_firm"] = wilson(sum(d["neg_control"]["class"] == "firm_positive" for d in run), len(run))
    out["neg_control_abs_delta_p95"] = tr(float(np.percentile([abs(d["neg_control"]["delta"]) for d in run], 95)))
    pos = [d for d in run if "delta" in d["pos_control"]]
    out["pos_control_firm"] = wilson(sum(d["pos_control"]["class"] == "firm_positive" for d in pos), len(pos))
    out["pos_control_point"] = wilson(sum(d["pos_control"]["class"] in ("firm_positive", "point_positive") for d in pos), len(pos))
    sc = [abs(d["sample_change"]) for d in run]
    out["sample_change_abs"] = {"median": tr(float(np.median(sc))), "max": tr(float(max(sc))),
                                "n_gt_log110": sum(v > C.LOG110 for v in sc)}
    d1 = [d for d in unfl if d["dag_family"] == "D1"]
    out["TyG_D1_firm_A"] = [sum(fp(d, "E3a_A") for d in d1), len(d1)]
    out["weight_kinds"] = {k: sum(d["weight_kind"] == k for d in run) for k in ("fasting", "mec", "diet_day1")}
    out["boot_failed_fits_total"] = sum(d["boot_failed_fits"] for d in run)
    out["N_h_star"] = {h: tr(v) for h, v in NHS.items()}
    json.dump(out, open(os.path.join(C.X, "m3_summary_i139.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


import csv
NEAR = {r["rec_id"] for r in csv.DictReader(open(os.path.join(C.X, "m1_results_i139.tsv"), encoding="utf-8"), delimiter="\t")
        if r["rank"] in ("strict_pass", "near_pass")}

if __name__ == "__main__":
    main()

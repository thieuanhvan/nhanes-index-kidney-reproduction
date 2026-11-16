"""I139 M3 two-order decomposition (prereg_i139_v7_en.md Section 6, M3, "Two-order decomposition"; descriptive only,
listed as not yet run in extraction/m2_m3_rules_i139.md rule 12). Written 15/11/2026 after all M3 results were seen.

From the paper's set P to the DAG A set by two routes, on the M3 common sample (rule 5: complete on y, the paper's terms,
the variables missing under DAG A and under DAG B, correct weight > 0):
  route 1 (add, then drop):  P -> P + addA -> P + addA - K
  route 2 (drop, then add):  P -> P - K    -> P - K + addA
addA = the variables missing under DAG A (same coding as m3_run_i139); K = the paper's terms mapped to the K-C, K-M,
K-O, K-X concepts of labels_i139.tsv (same term mapping as the E3b model of m3_desc_i139). Both routes end at the
same model (the DAG A set), so total = add step + drop step on either route. Reported: the add and drop steps on each
route and their mean over the two routes. Interpretation of "the DAG A set" (prereg silent): the paper's own terms minus
K plus the missing DAG-A variables, so precision and contested terms the paper adjusted for are kept. Uncertainty:
paired Rao-Wu bootstrap, B = 100 (descriptive; prereg sets no B for this analysis), seed 20261114-i139-twoorder:<rec_id>.
Usage: python -I tools/m3_twoorder_i139.py <rec_id> [...]  -> extraction/m3_twoorder/<rec_id>.json
Summary: python -I tools/m3_twoorder_i139.py --summary   -> extraction/m3_twoorder_summary_i139.json
"""
import os, sys, re, json, glob, math, time, traceback
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C
from m3_run_i139 import cls, concepts
from m3_desc_i139 import term_concept

B = int(os.environ.get("I139_TWOORDER_B", "100"))
OUT = os.path.join(C.X, "m3_twoorder")
os.makedirs(OUT, exist_ok=True)
tr = lambda v: math.floor(v * 1e4) / 1e4 if v is not None and np.isfinite(v) else v


def run(rid, T, L):
    t0 = time.time()
    lab = L[rid]
    b = C.load_builder(rid)
    if b.get("not_testable"):
        return {"rec_id": rid, "status": "not_testable"}
    df, y = b["df"], b["y"]
    dom = b["domain"].fillna(False).astype(bool)
    terms, tgt = list(b["terms"]), b["target"]
    w, kind = C.correct_weight(rid, b)
    df["_wc"] = w
    addA, addB = [], []
    for c in concepts(lab["missing_A"]):
        addA += [t for t in C.add_concept(df, c) if t not in terms and t not in addA]
    for c in concepts(lab["missing_B"]):
        addB += [t for t in C.add_concept(df, c) if t not in terms and t not in addB]
    stem = re.split(r"_q\d|__|_t\d|_seg", tgt)[0]
    protected = {t for t in terms if t == tgt or t.startswith(stem + "_") or t.startswith(stem + "__") or t == stem}
    kset = set()
    for col in ("adjusts_component", "adjusts_mediator", "conditions_on_outcome", "conditions_on_collider"):
        kset |= set(concepts(lab[col]))
    tmap = {t: term_concept(t) for t in terms if t not in protected}
    kterms = [t for t, c in tmap.items() if c in kset]
    allv = [y] + terms + addA + addB
    common = dom & df[allv].notna().all(axis=1) & (w > 0)
    tz, sd = C.standardise(df, tgt, common, w)
    P = [tz if t == tgt else t for t in terms]
    models = {"P": P, "A": P + addA, "D": [t for t in P if t not in kterms]}
    models["DA"] = models["D"] + addA
    mm = common.to_numpy()
    ww = w[mm].to_numpy(float)
    yy = df.loc[mm, y].to_numpy(float)
    log, prep = [], {}
    for k, ts in models.items():
        keep = [t for t in ts if t == tz or df.loc[common, t].nunique() > 1]
        if len(keep) < len(ts):
            log.append(f"model {k}: constant terms dropped {sorted(set(ts) - set(keep))}")
        Xm = np.column_stack([np.ones(mm.sum()), df.loc[mm, keep].to_numpy(float)])
        b0, _ = C.irls(Xm, yy, ww)
        prep[k] = (Xm, b0, 1 + keep.index(tz))
    point = {k: float(v[1][v[2]]) for k, v in prep.items()}
    keys, Mx = C.raowu_multipliers(df[df.SDMVSTRA.notna()], B, f"20261114-i139-twoorder:{rid}")
    kpos = {k: i for i, k in enumerate(keys)}
    rk = df.SDMVSTRA.fillna(-1).to_numpy(np.int64)[mm] * 1000 + df.SDMVPSU.fillna(-1).to_numpy(np.int64)[mm]
    ki = np.array([kpos[v] for v in rk])
    est = {k: np.full(B, np.nan) for k in prep}
    fails = 0
    for bb in range(B):
        wb = ww * Mx[bb, ki]
        for k, (Xm, b0, j) in prep.items():
            try:
                beta, _ = C.irls(Xm, yy, wb, beta0=b0, it=25, tol=1e-7)
                est[k][bb] = beta[j] if np.all(np.isfinite(beta)) else np.nan
            except np.linalg.LinAlgError:
                fails += 1

    def step(k1, k0):
        d = point[k1] - point[k0]
        v = est[k1] - est[k0]
        v = v[np.isfinite(v)]
        lo, hi = (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) > 20 else (np.nan, np.nan)
        return {"delta": d, "lo": lo, "hi": hi, "n_boot": int(len(v)), "class": cls(d, lo, hi)}

    def mean_step(a, b_):
        d = (point[a[0]] - point[a[1]] + point[b_[0]] - point[b_[1]]) / 2
        v = ((est[a[0]] - est[a[1]]) + (est[b_[0]] - est[b_[1]])) / 2
        v = v[np.isfinite(v)]
        lo, hi = (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) > 20 else (np.nan, np.nan)
        return {"delta": d, "lo": lo, "hi": hi, "class": cls(d, lo, hi)}

    return {"rec_id": rid, "status": "run", "weight_kind": kind, "n_common": int(mm.sum()), "B": B, "boot_failed_fits": fails,
            "added_A": addA, "k_terms_dropped": kterms, "log": log, "point_logor": point,
            "route1_add": step("A", "P"), "route1_drop": step("DA", "A"),
            "route2_drop": step("D", "P"), "route2_add": step("DA", "D"),
            "mean_add": mean_step(("A", "P"), ("DA", "D")), "mean_drop": mean_step(("DA", "A"), ("D", "P")),
            "total": step("DA", "P"), "seconds": int(time.time() - t0)}


def summary():
    R = [json.load(open(f, encoding="utf-8")) for f in sorted(glob.glob(os.path.join(OUT, "*.json")))]
    run_ = [d for d in R if d["status"] == "run"]
    out = {"n_files": len(R), "n_run": len(run_), "B": B, "errors": [d["rec_id"] for d in R if d["status"] == "error"]}
    for k in ("route1_add", "route1_drop", "route2_drop", "route2_add", "mean_add", "mean_drop", "total"):
        v = np.array([abs(d[k]["delta"]) for d in run_])
        sel = [d for d in run_ if abs(d[k]["delta"]) > 0]
        out[k] = {"n_nonzero": len(sel), "median_abs_all": tr(float(np.median(v))), "max_abs": tr(float(v.max())),
                  "median_abs_nonzero": tr(float(np.median([abs(d[k]["delta"]) for d in sel]))) if sel else None,
                  "n_gt_log110": int(sum(x > C.LOG110 for x in v)),
                  "n_firm": sum(d[k]["class"] == "firm_positive" for d in run_),
                  "n_point_only": sum(d[k]["class"] == "point_positive" for d in run_)}
    both = [d for d in run_ if d["added_A"] and d["k_terms_dropped"]]
    out["n_with_both_steps"] = len(both)
    out["order_dependence_add"] = {"median_abs": tr(float(np.median([abs(d["route1_add"]["delta"] - d["route2_add"]["delta"]) for d in both]))) if both else None,
                                   "max_abs": tr(float(max(abs(d["route1_add"]["delta"] - d["route2_add"]["delta"]) for d in both))) if both else None}
    out["share_total_from_drop_median"] = tr(float(np.median([abs(d["mean_drop"]["delta"]) / (abs(d["mean_drop"]["delta"]) + abs(d["mean_add"]["delta"]))
                                                               for d in run_ if abs(d["mean_drop"]["delta"]) + abs(d["mean_add"]["delta"]) > 0])))
    out["per_paper"] = {d["rec_id"]: {k: tr(d[k]["delta"]) for k in ("route1_add", "route1_drop", "route2_drop", "route2_add", "mean_add", "mean_drop", "total")}
                        for d in run_}
    json.dump(out, open(os.path.join(C.X, "m3_twoorder_summary_i139.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "per_paper"}, ensure_ascii=False, indent=1))


def main():
    if sys.argv[1:] == ["--summary"]:
        return summary()
    T, L = C.spec_table(), C.labels()
    for rid in sys.argv[1:]:
        try:
            res = run(rid, T, L)
        except Exception as e:
            traceback.print_exc()
            res = {"rec_id": rid, "status": "error", "error": f"{type(e).__name__}: {e}"[:400]}
        json.dump(res, open(os.path.join(OUT, rid + ".json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
        print(rid, res["status"], "add1", round(res.get("route1_add", {}).get("delta", np.nan), 4), "drop2",
              round(res.get("route2_drop", {}).get("delta", np.nan), 4), "total", round(res.get("total", {}).get("delta", np.nan), 4),
              res.get("seconds"), "s", flush=True)


if __name__ == "__main__":
    main()

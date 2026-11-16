"""I139 summary of the M3 descriptive outputs (tools/m3_desc_i139.py): E3b, PR / marginal-OR class changes, RoE,
conventional set and drug additions, as reported in verdict_vi.md 9.11 and the P6 main 5.6-5.7. Written 15/11/2026 for
the numexpr rerun (extraction/m234_correction_20261115.md); reads JSON only, refits nothing.
Usage: python -I tools/m3_desc_summary_i139.py <m3_desc dir> [<m3 dir>]
"""
import os, sys, json, glob, math
import numpy as np

LOG110 = math.log(1.10)
tr = lambda v: math.floor(v * 1e4) / 1e4


def main():
    dd = sys.argv[1]
    m3d = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(dd.rstrip("/\\")), "m3")
    R = [json.load(open(f, encoding="utf-8")) for f in sorted(glob.glob(os.path.join(dd, "*.json")))]
    R = [d for d in R if d.get("status") == "run"]
    M = {d["rec_id"]: d for d in (json.load(open(f, encoding="utf-8")) for f in glob.glob(os.path.join(m3d, "*.json")))}
    out = {"n_run": len(R)}
    e = [d for d in R if "E3b" in d]
    dl = np.array([d["E3b"]["logor"]["delta"] for d in e])
    away = sum(abs(d["point"]["E3b"]["logor"]) > abs(d["point"]["P"]["logor"]) for d in e)
    cl = [d["E3b"]["logor"]["class"] for d in e]
    out["E3b"] = {"n_with_terms": len(e), "away_from_1": int(away), "median_abs": tr(float(np.median(abs(dl)))),
                  "max_abs": tr(float(abs(dl).max())), "n_gt_log110": int((abs(dl) > LOG110).sum()),
                  "firm": cl.count("firm_positive"), "point": cl.count("point_positive")}
    for sc in ("PR", "marginal"):
        ch = []
        for d in R:
            for dag in ("A", "B"):
                k = f"E3a_{dag}_{sc}"
                if k in d and d["rec_id"] in M and d[k]["class"] != M[d["rec_id"]][f"E3a_{dag}"]["class"]:
                    ch.append(f'{d["rec_id"]} {dag}: {M[d["rec_id"]][f"E3a_{dag}"]["class"]} -> {d[k]["class"]}')
        out[f"class_changes_{sc}"] = ch
        out[f"firm_{sc}"] = sorted(d["rec_id"] for d in R for dag in ("A", "B")
                                   if f"E3a_{dag}_{sc}" in d and d[f"E3a_{dag}_{sc}"]["class"] == "firm_positive")
    roe = np.array([d["RoE"]["roe"] for d in R])
    out["RoE"] = {"median": tr(float(np.median(roe))), "q1": tr(float(np.percentile(roe, 25, method="hazen"))),
                  "q3": tr(float(np.percentile(roe, 75, method="hazen"))),   # quartile rule of the recorded summary
                  "max": tr(float(roe.max())), "n_cross_1": int(sum(0 < d["RoE"]["share_or_gt1"] < 1 for d in R))}
    for k in ("conv", "P_dmrx", "P_bprx"):
        v = np.array([abs(d[f"shift_{k}"]) for d in R])
        out[f"shift_{k}"] = {"median_abs": tr(float(np.median(v))), "max_abs": tr(float(v.max())), "n_gt_log110": int((v > LOG110).sum())}
    out["n_common"] = {d["rec_id"]: d["n_common"] for d in R}
    print(json.dumps({k: v for k, v in out.items() if k != "n_common"}, indent=1))
    return out


if __name__ == "__main__":
    main()

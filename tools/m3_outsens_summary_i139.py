"""I139 summary of the outcome-definition sensitivity runs (tools/m3_outsens_i139.py) against the main M3 run:
paper-model |change in log OR|, change in Delta_A, and class changes (A and B). Written 15/11/2026 for the numexpr
rerun (extraction/m234_correction_20261115.md); reads JSON only, refits nothing.
Usage: python -I tools/m3_outsens_summary_i139.py <m3_outsens dir> <m3 dir>
"""
import os, sys, json, glob, math
import numpy as np

LOG110 = math.log(1.10)
tr = lambda v: math.floor(v * 1e4) / 1e4


def main():
    od, m3d = sys.argv[1], sys.argv[2]
    M = {d["rec_id"]: d for d in (json.load(open(f, encoding="utf-8")) for f in glob.glob(os.path.join(m3d, "*.json")))}
    for var in sorted(os.listdir(od)):
        R = [json.load(open(f, encoding="utf-8")) for f in sorted(glob.glob(os.path.join(od, var, "*.json")))]
        R = [d for d in R if d.get("status") == "run" and d["rec_id"] in M]
        if not R:
            print(var, "no run")
            continue
        dP = np.array([abs(d["full"]["P"]["logor"] - M[d["rec_id"]]["full"]["P"]["logor"]) for d in R])
        dA = np.array([abs(d["E3a_A"]["delta"] - M[d["rec_id"]]["E3a_A"]["delta"]) for d in R
                       if "delta" in d["E3a_A"] and "delta" in M[d["rec_id"]]["E3a_A"]])
        ch = [f'{d["rec_id"]} {g}: {M[d["rec_id"]][f"E3a_{g}"]["class"]} -> {d[f"E3a_{g}"]["class"]}'
              for d in R for g in ("A", "B") if d[f"E3a_{g}"]["class"] != M[d["rec_id"]][f"E3a_{g}"]["class"]]
        firm = sorted(d["rec_id"] for d in R for g in ("A", "B") if d[f"E3a_{g}"]["class"] == "firm_positive")
        print(var, "papers", len(R), "| paper model |dlogOR| median", tr(float(np.median(dP))), "max", tr(float(dP.max())),
              "n>log1.10", int((dP > LOG110).sum()), "| change in Delta_A median", tr(float(np.median(dA))) if len(dA) else "na",
              "max", tr(float(dA.max())) if len(dA) else "na", "| class changes", ch, "| firm", firm)


if __name__ == "__main__":
    main()

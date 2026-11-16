"""I139 simulated positive control for the E3a comparison (rule 16, extraction/m2_m3_rules_i139.md, addition 3).
Usage: python -I tools/m3_simpc_i139.py <rec_id> [...]   -> extraction/m3_simpc/<rec_id>.json
"""
import os, sys, json, hashlib, time, traceback
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C
from m3_run_i139 import cls, concepts

B = 200
LEVELS = [1.05, 1.10, 1.20]
G = float(np.log(2.0))
OUT = os.path.join(C.X, "m3_simpc")
os.makedirs(OUT, exist_ok=True)


def rng_for(s):
    return np.random.default_rng(int(hashlib.sha256(s.encode()).hexdigest()[:15], 16))


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
    add = []
    for c in concepts(lab["missing_A"]) + concepts(lab["missing_B"]):
        add += [t for t in C.add_concept(df, c) if t not in terms and t not in add]
    common = dom & df[[y] + terms + add].notna().all(axis=1) & (w > 0)
    tz, _ = C.standardise(df, tgt, common, w)
    ts = [tz if t == tgt else t for t in terms]
    ts = [t for t in ts if t == tz or df.loc[common, t].nunique() > 1]
    mm = common.to_numpy()
    X0 = np.column_stack([np.ones(mm.sum()), df.loc[mm, ts].to_numpy(float)])
    yr = df.loc[mm, y].to_numpy(float)
    ww = w[mm].to_numpy(float)
    j = 1 + ts.index(tz)
    beta, _ = C.irls(X0, yr, ww)
    eta = X0 @ beta
    zx = X0[:, j]
    zx = (zx - np.average(zx, weights=ww)) / np.sqrt(np.average((zx - np.average(zx, weights=ww)) ** 2, weights=ww))
    keys, Mx = C.raowu_multipliers(df[df.SDMVSTRA.notna()], B, f"20261114-i139-simpc-boot:{rid}")
    kpos = {k: i for i, k in enumerate(keys)}
    rowkey = df.SDMVSTRA.fillna(-1).to_numpy(np.int64)[mm] * 1000 + df.SDMVPSU.fillna(-1).to_numpy(np.int64)[mm]
    ki = np.array([kpos[v] for v in rowkey])
    res = {"rec_id": rid, "status": "run", "n": int(mm.sum()), "levels": {}}
    for Lv in LEVELS:
        g = rng_for(f"20261114-i139-simpc:{rid}:{Lv}")
        e = g.standard_normal(mm.sum())
        e = e - np.average(e * zx, weights=ww) / np.average(zx * zx, weights=ww) * zx   # weighted-orthogonal to x
        e = (e - np.average(e, weights=ww)) / np.sqrt(np.average((e - np.average(e, weights=ww)) ** 2, weights=ww))
        u = g.random(mm.sum())

        def sim(r):
            Cv = r * zx + np.sqrt(1 - r * r) * e
            ys = (u < 1 / (1 + np.exp(-(eta + G * Cv)))).astype(float)
            bP, _ = C.irls(X0, ys, ww)
            XC = np.column_stack([X0, Cv])
            bC, _ = C.irls(XC, ys, ww)
            return bP[j] - bC[j], ys, XC, bP, bC
        lo, hi = -0.95, 0.95
        bias_hi, bias_lo = sim(hi)[0], sim(lo)[0]
        target = np.log(Lv)
        if not (bias_lo <= target <= bias_hi):
            res["levels"][str(Lv)] = {"class": "unreachable", "bias_range": [float(bias_lo), float(bias_hi)]}
            continue
        for _ in range(20):
            mid = (lo + hi) / 2
            if sim(mid)[0] < target:
                lo = mid
            else:
                hi = mid
        r = (lo + hi) / 2
        bias, ys, XC, bP, bC = sim(r)
        d_full = bC[j] - bP[j]
        v = np.full(B, np.nan)
        for bb in range(B):
            wb = ww * Mx[bb, ki]
            try:
                p1, _ = C.irls(X0, ys, wb, beta0=bP, it=25, tol=1e-7)
                p2, _ = C.irls(XC, ys, wb, beta0=bC, it=25, tol=1e-7)
                v[bb] = p2[j] - p1[j]
            except np.linalg.LinAlgError:
                pass
        v = v[np.isfinite(v)]
        lo_, hi_ = np.percentile(v, 2.5), np.percentile(v, 97.5)
        res["levels"][str(Lv)] = {"r": float(r), "delta": float(d_full), "lo": float(lo_), "hi": float(hi_),
                                  "n_boot": int(len(v)), "class": cls(d_full, lo_, hi_)}
    res["seconds"] = int(time.time() - t0)
    return res


def main():
    T, L = C.spec_table(), C.labels()
    for rid in sys.argv[1:]:
        try:
            res = run(rid, T, L)
        except Exception as e:
            traceback.print_exc()
            res = {"rec_id": rid, "status": "error", "error": f"{type(e).__name__}: {e}"[:300]}
        json.dump(res, open(os.path.join(OUT, rid + ".json"), "w", encoding="utf-8"), indent=1, default=float)
        print(rid, res["status"], {k: v.get("class") for k, v in res.get("levels", {}).items()}, res.get("seconds"), flush=True)


if __name__ == "__main__":
    main()

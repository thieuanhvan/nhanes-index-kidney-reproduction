"""I139: which single-file numexpr fault reproduces the n of the M2 rows recorded at 08:37-08:39 (before the fix)?
For one paper: load the builder once (fixed code) to fill the reader caches, then for every cached XPT file with more
than 1e6 numeric cells (the size at which pandas hands DataFrame.mask to numexpr) replace the cached frame by the
all-zero frame that the fault produces (every numeric value 0, SEQN 0), rebuild the paper and record n of the M2
regimes: correct weight (i, iv), wrong-type weight (ii), unweighted (iii, v). Restore the file and go on.
Usage: python -I tools/m2_fileloss_sim_i139.py <rec_id> <n_i> <n_ii> <n_iii> [...]  (recorded n values to match)
Output: one line per file, '<-- MATCH' when all three n equal the recorded ones.
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import m23_common_i139 as C
lib, nc, nx = C.lib, C.nc, C.nx


def ns(rid):
    lib._tab.clear()
    b = C.load_builder(rid)
    df = b["df"]
    dom = b["domain"].fillna(False).astype(bool)
    cc = dom & df[[b["y"]] + b["terms"]].notna().all(axis=1)
    wc, kind = C.correct_weight(rid, b)
    ww = C.wrong_weight(df, kind)
    return int((cc & (wc > 0)).sum()), int((cc & (ww > 0)).sum()), int(cc.sum()), kind


def zeroed(x, indexed):
    if indexed:   # nhanes_extra cache: drop_duplicates('SEQN').set_index('SEQN') after the fault
        z = x.reset_index()
        num = z.select_dtypes("number").columns
        z[num] = 0.0
        z["SEQN"] = 0
        return z.drop_duplicates("SEQN").set_index("SEQN")
    z = x.copy()
    num = z.select_dtypes("number").columns
    z[num] = 0.0
    z["SEQN"] = z["SEQN"].astype(np.int64)
    return z


def main():
    rid = sys.argv[1]
    target = tuple(int(v) for v in sys.argv[2:5])
    base = ns(rid)
    print(rid, "fixed code n (i, ii, iii)", base[:3], "kind", base[3], "recorded", target, flush=True)
    caches = [("nhanes_core.read", nc._cache, False), ("nhanes_extra.read_stem", nx._cache, True), ("m1_lib._read", lib._raw, False)]
    for name, cache, indexed in caches:
        for p in list(cache):
            x = cache[p]
            cells = x.shape[0] * (x.select_dtypes("number").shape[1] + (1 if indexed else 0))
            if cells <= 1_000_000:
                continue
            cache[p] = zeroed(x, indexed)
            try:
                r = ns(rid)
            finally:
                cache[p] = x
            tag = "  <-- MATCH" if r[:3] == target else ""
            print(name, os.path.relpath(p, nc.ROOT), cells, "n", r[:3], r[3], tag, flush=True)


if __name__ == "__main__":
    main()

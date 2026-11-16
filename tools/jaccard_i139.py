"""I139 overlap of people between papers (prereg Section 4): Jaccard index on the SEQN sets of the reconstructed
analytic samples (builder domain, complete case on y and the paper's terms, correct weight > 0), for each pair of
papers with overlapping cycles. SEQN is unique across NHANES cycles.
Usage: python -I tools/jaccard_i139.py  -> extraction/jaccard_i139.tsv and summary printed
"""
import os, sys, csv, json, itertools
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C


def main():
    T = C.spec_table()
    sets, cyc = {}, {}
    cache = os.path.join(C.X, "jaccard_sets_i139.json")
    if os.path.exists(cache):
        d = json.load(open(cache))
        sets = {k: set(v) for k, v in d["sets"].items()}; cyc = {k: set(v) for k, v in d["cyc"].items()}
    for rid in T:
        if rid in sets:
            continue
        b = C.load_builder(rid)
        if b.get("not_testable"):
            continue
        df = b["df"]
        w, _ = C.correct_weight(rid, b)
        m = b["domain"].fillna(False).astype(bool) & df[[b["y"]] + b["terms"]].notna().all(axis=1) & (w > 0)
        sets[rid] = set(int(s) for s in df.loc[m, "SEQN"])
        cyc[rid] = set(df.loc[m, "cycle"].unique())
        print(rid, len(sets[rid]), flush=True)
        json.dump({"sets": {k: sorted(v) for k, v in sets.items()}, "cyc": {k: sorted(v) for k, v in cyc.items()}}, open(cache, "w"))
    rows = []
    for a, b in itertools.combinations(sorted(sets), 2):
        if not (cyc[a] & cyc[b]):
            continue
        i = len(sets[a] & sets[b]); u = len(sets[a] | sets[b])
        rows.append({"a": a, "b": b, "family_a": T[a]["family"], "family_b": T[b]["family"], "jaccard": int(i / u * 10000) / 10000,
                     "overlap_n": i, "share_of_smaller": int(i / min(len(sets[a]), len(sets[b])) * 10000) / 10000})
    with open(os.path.join(C.X, "jaccard_i139.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    j = np.array([r["jaccard"] for r in rows]); s = np.array([r["share_of_smaller"] for r in rows])
    print("pairs with overlapping cycles", len(rows), "median Jaccard", np.median(j), "p90", np.percentile(j, 90), "max", j.max(),
          "pairs Jaccard > 0.5", int((j > 0.5).sum()), "pairs share_of_smaller > 0.8", int((s > 0.8).sum()))


if __name__ == "__main__":
    main()

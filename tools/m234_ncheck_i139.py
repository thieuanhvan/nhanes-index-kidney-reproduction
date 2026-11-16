"""I139 sample check after the numexpr fix (extraction/m234_correction_20261115.md).
For every testable paper: load the builder with the fixed loaders, take the correct weight (rule 1) and record
n of the paper sample (builder domain, complete on y and the paper's terms, correct weight > 0; the M3 n_paper_sample
and the Jaccard set), n per cycle, the cycles of the builder frame with no row in the sample, and a sha256 of the
sorted SEQN set. Run in a separate process from M2/M3/M4/Jaccard; agreement between independent loads, and no cycle of
the frame missing from the sample, show that no cycle was lost.
Usage: python -I tools/m234_ncheck_i139.py <out.tsv>
"""
import os, sys, csv, json, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import m23_common_i139 as C


def main():
    assert pd.get_option("compute.use_numexpr") is False, "numexpr still on"
    T = C.spec_table()
    rows = []
    for rid in sorted(T):
        b = C.load_builder(rid)
        if b.get("not_testable"):
            continue
        df = b["df"]
        w, kind = C.correct_weight(rid, b)
        m = b["domain"].fillna(False).astype(bool) & df[[b["y"]] + b["terms"]].notna().all(axis=1) & (w > 0)
        per = df.loc[m, "cycle"].value_counts().to_dict()
        frame_cyc = sorted(set(df.cycle))
        empty = [c for c in frame_cyc if per.get(c, 0) == 0]
        seqn = sorted(int(s) for s in df.loc[m, "SEQN"])
        rows.append({"rec_id": rid, "weight_kind": kind, "n_paper_sample": int(m.sum()),
                     "n_per_cycle": json.dumps({c: int(per.get(c, 0)) for c in frame_cyc}),
                     "frame_cycles_empty": "|".join(empty),
                     "seqn_sha256": hashlib.sha256(",".join(map(str, seqn)).encode()).hexdigest()})
        print(rid, kind, int(m.sum()), "empty:", empty, flush=True)
    with open(sys.argv[1], "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()

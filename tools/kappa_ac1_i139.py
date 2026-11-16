"""I139 prereg_i139_v7_en.md Section 3.2 step 1: per-field raw agreement, Cohen's kappa and Gwet's AC1 between the two
independent extraction passes (extraction/pass_A, extraction/pass_B), main extraction of 55 papers. Written 15/11/2026.
Values are normalised exactly as in tools/compare_extraction_i139.py (same fields, same format-only equivalences), so the
raw agreement reproduces the 845/1045 printed by that script. Each distinct normalised value is a nominal category.
Cohen's kappa: (po - pe) / (1 - pe), pe = sum_c pA(c) pB(c) (each pass's own marginals).
Gwet's AC1: (po - pe) / (1 - pe), pe = sum_c pi_c (1 - pi_c) / (q - 1), pi_c = mean of the two passes' marginals, q = number
of categories observed (same formula as compare_extraction_i139.py). kappa is undefined (reported NA) when pe = 1.
Numbers truncated to 4 decimals.
Usage: python -I tools/kappa_ac1_i139.py   -> extraction/kappa_ac1_i139.tsv (and printed)
"""
import os, sys, json, csv, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_extraction_i139 import KEY, FREE, get, norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X = os.path.join(ROOT, "extraction")
tr = lambda v: "NA" if v is None else str(math.floor(v * 1e4) / 1e4)


def stats(pairs):
    n = len(pairs)
    po = sum(a == b for a, b in pairs) / n
    cats = sorted({x for p in pairs for x in p})
    pA = {c: sum(a == c for a, _ in pairs) / n for c in cats}
    pB = {c: sum(b == c for _, b in pairs) / n for c in cats}
    pe_k = sum(pA[c] * pB[c] for c in cats)
    kappa = (po - pe_k) / (1 - pe_k) if pe_k < 1 else None
    if len(cats) > 1:
        pis = [(pA[c] + pB[c]) / 2 for c in cats]
        pe_g = sum(p * (1 - p) for p in pis) / (len(cats) - 1)
    else:
        pe_g = 0.0
    ac1 = (po - pe_g) / (1 - pe_g) if pe_g < 1 else 1.0
    return po, kappa, ac1, len(cats)


def main():
    da, db = os.path.join(X, "pass_A"), os.path.join(X, "pass_B")
    ids = sorted({f[:-5] for f in os.listdir(da) if f.endswith(".json")} & {f[:-5] for f in os.listdir(db) if f.endswith(".json")})
    A = {i: json.load(open(os.path.join(da, i + ".json"), encoding="utf-8")) for i in ids}
    Bp = {i: json.load(open(os.path.join(db, i + ".json"), encoding="utf-8")) for i in ids}
    rows, ag, tot = [], 0, 0
    for f in KEY:
        pairs = [(norm(get(A[i], f), f), norm(get(Bp[i], f), f)) for i in ids]
        po, k, ac1, q = stats(pairs)
        a = sum(x == y for x, y in pairs)
        ag += a; tot += len(pairs)
        rows.append([f, "free text" if f in FREE else "", a, len(pairs), tr(po), tr(k), tr(ac1), q])
    with open(os.path.join(X, "kappa_ac1_i139.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["field", "note", "agree", "n", "raw_agreement", "cohen_kappa", "gwet_ac1", "n_categories"])
        w.writerows(rows)
    print(f"papers {len(ids)}; overall raw agreement {ag}/{tot}")
    for r in rows:
        print("\t".join(map(str, r)))


if __name__ == "__main__":
    main()

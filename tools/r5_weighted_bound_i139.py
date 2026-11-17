"""I139 external review round 5 (18/11/2026): upper bound for the N_h*-weighted primary E3a when the stratified
bootstrap is degenerate (every stratum has 0 firm positives). Post hoc and descriptive, no refit.

Effective sample size of the stratified estimator (Korn and Graubard): n_eff = 1 / sum_h(W_h^2 / n_h), with the
weights W_h of the registered analysis (m3_summary_i139.NHS). With 0 events the upper 95 % bound is the Wilson bound
z^2 / (n_eff + z^2) and the Clopper-Pearson bound 1 - 0.025^(1 / n_eff).
Usage: python -I tools/r5_weighted_bound_i139.py   -> extraction/r5_weighted_bound_i139.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_summary_i139 as S

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X = os.path.join(ROOT, "extraction")


def main():
    summ = json.load(open(os.path.join(X, "m3_summary_i139.json"), encoding="utf-8"))
    out = {}
    for key in ("E3a_A_primary", "E3a_B_primary"):
        by = summ[key]["by_stratum"]
        k = sum(v[0] for v in by.values())
        tot = sum(S.NHS[h] for h in by)
        W = {h: S.NHS[h] / tot for h in by}
        neff = 1 / sum(W[h] ** 2 / by[h][1] for h in by)
        rec = {"events": k, "n": sum(v[1] for v in by.values()), "n_eff": S.tr(neff),
               "weights": {h: S.tr(W[h]) for h in sorted(W)}}
        if k == 0:
            rec["wilson_upper"] = S.tr(S.Z ** 2 / (neff + S.Z ** 2))
            rec["clopper_pearson_upper"] = S.tr(1 - 0.025 ** (1 / neff))
        out[key] = rec
    json.dump(out, open(os.path.join(X, "r5_weighted_bound_i139.json"), "w", encoding="utf-8"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

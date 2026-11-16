"""I139 M2 (survey design regimes) on M1 strict/near-pass papers, per extraction/m2_m3_rules_i139.md rules 13-15.
Usage: python -I tools/m2_run_i139.py   (reads extraction/m1_results_i139.tsv; writes extraction/m2_results_i139.tsv)
"""
import os, sys, csv, math
import numpy as np, pandas as pd
from scipy import stats
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C

OUTP = os.path.join(C.X, "m2_results_i139.tsv")
T = lambda v: math.floor(v * 1e4) / 1e4 if np.isfinite(v) else v


def model_se_fit(df, y, terms, w, dom):
    cc = dom & df[[y] + terms].notna().all(axis=1) & (w > 0)
    Xm = np.column_stack([np.ones(cc.sum()), df.loc[cc, terms].to_numpy(float)])
    yy = df.loc[cc, y].to_numpy(float)
    ww = w[cc].to_numpy(float)
    beta, H = C.irls(Xm, yy, ww)
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    return beta, se, int(cc.sum())


def main():
    m1 = [r for r in csv.DictReader(open(os.path.join(C.X, "m1_results_i139.tsv"), encoding="utf-8"), delimiter="\t")
          if r["rank"] in ("strict_pass", "near_pass")]
    S = C.spec_table()
    rows = []
    for r in m1:
        rid = r["rec_id"]
        sp = S[rid]
        b = C.load_builder(rid)
        df, y, terms, tgt = b["df"], b["y"], list(b["terms"]), b["target"]
        dom = b["domain"].fillna(False).astype(bool)
        wc, kind = C.correct_weight(rid, b)
        df["_wc"] = wc
        df["_ww"] = C.wrong_weight(df, kind)
        df["_w1"] = 1.0
        j = 1 + terms.index(tgt)
        reg = {}
        for name, wcol, se_type in [("i_correct_design", "_wc", "design"), ("ii_wrongtype_design", "_ww", "design"),
                                    ("iii_unweighted_model", "_w1", "model"), ("iv_weighted_model", "_wc", "model"),
                                    ("v_unweighted_design", "_w1", "design")]:
            if se_type == "design":
                f = C.nc.svy_logit(df, y, terms, wcol, domain=dom)
                t = f["table"].loc[tgt]
                lo, hi, orr, n = float(t.lo), float(t.hi), float(t.OR), f["n"]
            else:
                beta, se, n = model_se_fit(df, y, terms, df[wcol], dom)
                orr = math.exp(beta[j]); lo = math.exp(beta[j] - 1.959963984540054 * se[j]); hi = math.exp(beta[j] + 1.959963984540054 * se[j])
            reg[name] = (orr, lo, hi, n)
        o, l, h, nr = float(sp["estimate"]), float(sp["ci_low"]), float(sp["ci_high"]), float(sp["n_analytic"])
        wr = math.log(h) - math.log(l)
        comp = {}
        for k, (orr, lo, hi, n) in reg.items():
            width = (math.log(hi) - math.log(lo)) * math.sqrt(n / nr)
            comp[k] = (abs(math.log(orr) - math.log(o)), abs(math.log(width) - math.log(wr)) if wr > 0 else float("nan"))
        order_or = sorted(comp, key=lambda k: comp[k][0])
        order_w = sorted(comp, key=lambda k: comp[k][1])
        match = "cannot_be_distinguished"
        if wr <= 0:
            match = "cannot_be_distinguished_reported_ci_width_0"
        elif order_or[0] == order_w[0]:
            k0 = order_or[0]
            gap_or = comp[order_or[1]][0] - comp[k0][0]
            gap_w = comp[order_w[1]][1] - comp[k0][1]
            if gap_or >= C.LOG110 or gap_w >= C.LOG110:
                match = k0
        oi, loi, hii, ni = reg["i_correct_design"]
        ratio = wr / ((math.log(hii) - math.log(loi)) * math.sqrt(ni / nr)) if wr > 0 else float("nan")
        # OR per SD ratio (matched regime vs correct design): scale by weighted SD of a continuous target
        sdv = 1.0
        x = df.loc[dom & df[[y] + terms].notna().all(axis=1) & (wc > 0), tgt]
        if x.nunique() > 2:
            sdv = float(np.sqrt(np.cov(x, aweights=wc[x.index])))
        ormatch = reg[match][0] if match in reg else float("nan")
        or_sd_ratio = math.exp(sdv * (math.log(ormatch) - math.log(oi))) if match in reg else float("nan")
        row = {"rec_id": rid, "m1_rank": r["rank"], "weight_kind_correct": kind, "or_reported": o, "ci_reported": f"{l}-{h}",
               "n_reported": int(nr), "matched_regime": match, "width_ratio_E2": T(ratio), "E2_lt_080": (ratio < 0.80) if wr > 0 else "na",
               "or_per_sd_ratio_matched_vs_correct": T(or_sd_ratio)}
        for k, (orr, lo, hi, n) in reg.items():
            row[k] = f"OR {T(orr)} [{T(lo)}; {T(hi)}] n {n} dlogOR {T(comp[k][0])} dlogW {T(comp[k][1])}"
        rows.append(row)
        print(rid, r["rank"], kind, "match", match, "E2 ratio", T(ratio), flush=True)
    with open(OUTP, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()

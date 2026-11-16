"""KD3 / prereg_i139_v7_en.md section 3.2 step 9: variance implementation check.
(A) nhanes_core.svy_logit vs samplics.SurveyGLM (independent implementation) on 3 PILOT specifications
    (pilot data tyg_pooled_0518_v2.pkl, specs taken from pilot/i139_fix/r1_r2_specs.csv). Pass: every coefficient
    and SE within 1 % relative difference.
(B) reproduce a published NCHS estimate (Data Brief 516, 2021-2023 diabetes prevalence) with nhanes_core.
No frame paper is touched; no new index x kidney specification is fitted.
python -I kd3_variance_check.py  -> pilot/kd3_variance_check.json"""
import json, os, site, sys, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.append(site.getusersitepackages())          # samplics + polars installed with pip --user (-I hides user site)
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import nhanes_core as nc
import samplics
from samplics import SurveyGLM, ModelType

HERE = os.path.dirname(os.path.abspath(__file__)); PILOT = os.path.join(HERE, "..", "pilot")
out = {"samplics_version": samplics.__version__}

# ---------------- (A) pilot specifications ----------------
d = pd.read_pickle(os.path.join(PILOT, "tyg_pooled_0518_v2.pkl"))
d = d[(d.RIDAGEYR >= 20) & (d.RIDEXPRG != 1)].copy().reset_index(drop=True)          # as run_pilot_i139.py
for c in ["WTSAF2YR", "WTMEC2YR"]:   # XPT zero read as 5.4e-79 by pandas: exact 0 (changes N only, not estimates)
    d[c] = d[c].mask(d[c].abs() < 1e-70, 0.0)
d["ckd21"] = ((d.uacr >= 30) | (d.egfr < 60)).astype(float).where(d.uacr.notna() & d.egfr.notna())
d["albu_gt"] = (d.uacr > 30).astype(float).where(d.uacr.notna())
d["female"] = d.female.astype(float)
for c in ["race", "educ", "smoke"]:
    for lv in sorted(d[c].dropna().unique())[1:]:
        d[f"{c}_{lv}"] = (d[c] == lv).astype(float).where(d[c].notna())
COV = {"age": ["RIDAGEYR"], "sex": ["female"], "race": [c for c in d if c.startswith("race_")], "educ": [c for c in d if c.startswith("educ_")],
       "pir": ["INDFMPIR"], "smoke": [c for c in d if c.startswith("smoke_")], "htn": ["htn"], "dm": ["dm"], "bmi": ["BMXBMI"], "scr": ["LBXSCR"],
       "egfr09": ["egfr09"], "tc": ["LBXSCH"], "ua": ["LBXSUA"], "alb": ["LBXSAL"], "bun": ["LBXSBU"]}
ALL = [2005, 2007, 2009, 2011, 2013, 2015, 2017]
SPECS = [  # (label, paper, cycles, domain, y, exposure, covariates, weight)
    ("S1 B1 adults 2011-2018, albuminuria, fasting weight", "B1", [2011, 2013, 2015, 2017], "adult", "albu_gt", "tyg",
     ["age", "sex", "race", "smoke", "htn", "dm", "bmi", "scr", "egfr09"], "fast"),
    ("S2 B2 adults 2005-2018, CKD (2021 eq.), MEC weight, non-fasting index", "B2", ALL, "adult", "ckd21", "tyg_nonfast",
     ["age", "sex", "race", "pir", "smoke", "htn", "dm", "tc"], "mec"),
    ("S3 B4 diabetes domain 2007-2016, CKD (2021 eq.), fasting weight", "B4", [2007, 2009, 2011, 2013, 2015], "dm", "ckd21", "tyg",
     ["age", "sex", "race", "educ", "smoke", "htn", "pir", "bmi", "alb", "bun", "ua", "tc"], "fast"),
]
pilot_csv = pd.read_csv(os.path.join(PILOT, "i139_fix", "r1_r2_specs.csv"))
resA = []
for label, paper, cyc, dom, y, ex, covs, wm in SPECS:
    sub = d[d.cycle.isin(cyc)].reset_index(drop=True)
    x = [ex] + sum([COV[c] for c in covs], [])
    wcol = "WTSAF2YR" if wm == "fast" else "WTMEC2YR"
    sub["w"] = sub[wcol].fillna(0) / len(cyc)
    domain = ((sub.dm == 1) if dom == "dm" else pd.Series(True, index=sub.index)).to_numpy() & sub[ex].notna().to_numpy()
    r = nc.svy_logit(sub, y, x, "w", domain=domain)
    # samplics: all rows of the pooled design, rows outside the analysed set get weight 0 (domain estimation)
    cc = domain & sub[[y] + x].notna().all(axis=1).to_numpy() & (sub.w > 0).to_numpy()
    X = sub[x].to_numpy(float).copy(); X[~cc] = 0.0
    Y = np.where(cc, sub[y].fillna(0).to_numpy(float), 0.0); W = np.where(cc, sub.w.to_numpy(), 0.0)
    g = SurveyGLM(model=ModelType.LOGISTIC)
    g.estimate(y=Y, x=X, x_labels=x, samp_weight=W, stratum=sub.SDMVSTRA.to_numpy(int), psu=(sub.SDMVSTRA.to_numpy(int) * 10 + sub.SDMVPSU.to_numpy(int)), add_intercept=True)
    bs, ss = np.asarray(g.beta["point_est"]), np.asarray(g.beta["stderror"])
    rel_b = np.abs(r["coef"] - bs) / np.maximum(np.abs(bs), 1e-12); rel_s = np.abs(r["se"] - ss) / ss
    # pilot output for the same spec (OR of the exposure, design df = 4 sig. digits)
    pc = pilot_csv[(pilot_csv.paper == paper) & (pilot_csv.y == y) & (pilot_csv.wmode == wm) & (pilot_csv.expo == ex)].iloc[0]
    tab = r["table"].loc[ex]
    resA.append({"spec": label, "n": r["n"], "n_pilot": int(pc.n), "df": r["df"], "exposure": ex,
                 "coef_core": float(r["coef"][1]), "coef_samplics": float(bs[1]), "se_core": float(r["se"][1]), "se_samplics": float(ss[1]),
                 "OR_core": float(tab.OR), "lo_core": float(tab.lo), "hi_core": float(tab.hi), "OR_pilot": float(pc.OR), "lo_pilot": float(pc.lo), "hi_pilot": float(pc.hi),
                 "max_rel_diff_coef_all": float(rel_b.max()), "max_rel_diff_se_all": float(rel_s.max()), "n_params": len(bs),
                 "pass": bool(rel_b.max() <= 0.01 and rel_s.max() <= 0.01)})
    print(json.dumps(resA[-1]), flush=True)
out["A"] = resA

# ---------------- (B) NCHS Data Brief 516, August 2021-August 2023 ----------------
L = nc.build(["2021-2023"])
L["wf"] = nc.pooled_weight(L, "fasting")
base = (L.age >= 20) & (L.pregnant != 1) & (L.wf > 0)
fast8 = L.fast_hours.between(8, 24)
dx = L.DIQ010 == 1
und = (~dx) & L.DIQ010.isin([2, 3]) & ((L.glucose >= 126) | (L.hba1c >= 6.5))
valid = L.DIQ010.isin([1, 2, 3]) & L.glucose.notna() & L.hba1c.notna()
L["diag"] = dx.astype(float); L["undiag"] = und.astype(float); L["total"] = (dx | und).astype(float)
STD = {"20-39": 0.396579, "40-59": 0.371795, "60+": 0.231626}  # US 2000 standard, NCHS age groups 20-39/40-59/60+
L["ag"] = pd.cut(L.age, [19, 39, 59, 200], labels=list(STD))
resB = {}
for vname, dom in [("fasting 8-24 h, valid DIQ010+FPG+HbA1c", base & fast8 & valid), ("fasting weight > 0, valid DIQ010+FPG+HbA1c", base & valid)]:
    rr = {"n": int(dom.sum())}
    for k in ["diag", "undiag", "total"]:
        m = nc.svy_mean(L, k, "wf", domain=dom.to_numpy())
        aa = sum(STD[a] * nc.svy_mean(L, k, "wf", domain=(dom & (L.ag == a)).to_numpy())["est"] for a in STD)
        rr[k] = {"crude_pct": 100 * m["est"], "se_pct": 100 * m["se"], "df": m["df"], "age_adj_pct": 100 * aa}
    resB[vname] = rr
    print(vname, json.dumps(rr), flush=True)
out["B"] = resB
out["B_published"] = {"source": "NCHS Data Brief No. 516 (Nov 2024), https://www.cdc.gov/nchs/products/databriefs/db516.htm",
                      "crude": {"total": 15.8, "diag": 11.3, "undiag": 4.5}, "age_adj": {"total": 14.3, "diag": 10.1, "undiag": 4.2}}
json.dump(out, open(os.path.join(PILOT, "kd3_variance_check.json"), "w"), indent=1)

"""M1 builder, aca6f864fb (OBS and CKD; source withheld).
Spec: OBS quartiles, Q4 vs Q1, Model 3 (all covariates), CKD = eGFR < 60 or ACR > 30, NHANES 2005-2018, weighted (2-year weight / 7).
Sources (fulltext_txt/aca6f864fb.txt): exclusions section 2.1, outcome section 2.2 (CKD-EPI 2009 printed with race factor),
OBS section 2.3, covariates section 2.4, weighting and models section 2.5, Table 1 and Table 2.
OBS from lib.DIET["aca6f864fb"](df, reading=1) (extraction/diet_recon_notes_i139.md).
"""
import numpy as np
import pandas as pd

REC_ID = "aca6f864fb"
CYCLES = ["2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def labels(index, pairs):
    out = pd.Series(np.nan, index=index, dtype=object)
    for m, lab in reversed(pairs):
        out[m.fillna(False).to_numpy(bool)] = lab
    return out


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.DIET[REC_ID](df, reading=1)
    e = df.egfr_ckdepi2009
    df["y"] = ((e < 60) | (df.uacr > 30)).astype(float).where(e.notna() & df.uacr.notna())
    dom = (df.age >= 18) & df.scr.notna() & df.uacr.notna() & df.x.notna()
    qn, cuts = lib.quantile_dummies(df, "x", dom, k=4)
    # covariates (section 2.4)
    df["race4"] = df.race1.map({"mexican_american": "mex", "nh_white": "white", "nh_black": "black",
                                "other_hispanic": "other", "other_multi": "other"})
    mar = lib.raw(df, lib.stems("DEMO"), "DMDMARTL")
    df["marital3"] = mar.map({1: "married_partner", 6: "married_partner", 2: "dsw", 3: "dsw", 4: "dsw", 5: "never"})
    df["pir3"] = labels(df.index, [(df.pir < 1.3, "low"), (df.pir <= 3.5, "mid"), (df.pir > 3.5, "high")])
    sl = lib.raw(df, lib.stems("SLQ"), ["SLD010H", "SLD012"])
    df["sleep_h"] = sl.where(sl < 77)
    dpq = pd.DataFrame({i: lib.raw(df, lib.stems("DPQ"), f"DPQ0{i}0") for i in range(1, 10)})
    dpq = dpq.where(dpq <= 3)
    phq = dpq.sum(axis=1, min_count=9)
    df["depress"] = flag(phq >= 10, phq.notna())
    df["htn"] = df["hypertension_any_140"]
    pills = lib.raw(df, lib.stems("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, lib.stems("DIQ"), "DIQ050")
    ogtt = lib.raw(df, lib.stems("OGTT"), "LBXGLT")
    df["dm"] = flag((df.dx_diabetes == 1) | (pills == 1) | (ins == 1) | (df.hba1c >= 6.5) | (df.glucose >= 126)
                    | (df.glucose_serum >= 200) | (ogtt >= 200),
                    df.dx_diabetes.notna() | df.hba1c.notna() | df.glucose.notna())
    lipmed = lib.raw(df, lib.stems("BPQ"), "BPQ100D")
    low_hdl = pd.Series(np.where(df.female == 1, df.hdl <= 50, df.hdl <= 40), index=df.index) & df.hdl.notna()
    df["hlip"] = flag((lipmed == 1) | (df.tc >= 200) | (df.tg >= 150) | (df.ldl >= 130) | low_hdl,
                      df.tc.notna() & df.hdl.notna())
    df["dm_drug"] = flag((pills == 1) | (ins == 1), df.dx_diabetes.notna())
    bpmed = lib.raw(df, lib.stems("BPQ"), "BPQ050A")
    df["htn_drug"] = flag(bpmed == 1, df.dx_hypertension.notna())
    terms = list(qn)
    terms += ["age", "female"]
    terms += lib.dummies(df, "race4", "white")
    terms += lib.dummies(df, "marital3", "married_partner")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "pir3", "low")
    terms += ["dm", "hlip", "htn", "cvd_history", "depress", "sleep_h", "dm_drug", "htn_drug"]
    df["w"] = lib.pooled_weight(df, "mec")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x_q4", "weight": "w",
            "readings": {
                "ACR threshold (spec ambiguity)": "reading_1: ACR > 30 mg/g (strict)",
                "which survey weight (spec ambiguity)": "reading_1: 2-year weight / 7; kind = MEC (smallest subsample among the lib kinds that the analysis needs: urine and serum labs; dietary weight is reading_2)",
                "Model 2 vs Model 3 (spec ambiguity)": "reading_1: Model 3, all covariates",
                "eGFR equation": "CKD-EPI 2009 with race (printed formula with 1.159 black factor)",
                "OBS": "diet_indices_i139 reading 1 (mean of recalls, sex-specific tertiles among adults >= 18, score kept with >= 16 of 20 items)",
                "quartile cut": "unweighted OBS quartiles in the domain (age >= 18, Scr, UACR, OBS present); no weighted quantiles stated",
                "age": "continuous",
                "race": "4 groups as listed: NH white (reference), NH black, Mexican American, other (incl. other Hispanic)",
                "marital": "3 groups: divorced/separated/widowed, married/living with partner (reference), never married",
                "education": "3 groups DMDEDUC2 1-2 / 3 / 4-5",
                "PIR": "< 1.3, 1.3 to 3.5, > 3.5 (methods text)",
                "sleep duration": "hours, SLD010H (2005-2014) or SLD012 (2015-2018, weekday/workday), continuous",
                "depression": "PHQ-9 total >= 10, all 9 items required (codes 7, 9 missing)",
                "hypertension": "BPQ050A = 1, BPQ020 = 1, SBP >= 140 or DBP >= 90",
                "diabetes": "pills (DIQ070; DID070 in 2005-2008) or insulin (DIQ050), DIQ010 = 1, HbA1c >= 6.5, fasting glucose >= 126 mg/dL, random serum glucose (BIOPRO LBXSGL) >= 200 mg/dL (11.1 mmol/L), OGTT 2 h >= 200 mg/dL (2005-2016)",
                "hyperlipidemia": "BPQ100D = 1 (now taking cholesterol medicine), TC >= 200, TG >= 150, LDL >= 130, or HDL <= 50 women / <= 40 men; known when TC and HDL present",
                "CVD": "cvd_history (MCQ160B-F: CHF, CHD, angina, heart attack, stroke)",
                "hypoglycemic drug": "pills or insulin = 1, else 0 when DIQ010 answered",
                "antihypertensive drug": "BPQ050A = 1, else 0 when BPQ020 answered",
                "covariate exclusion (exclusion 4)": "complete case on model terms by the fit; quartiles are cut before it (on the domain), as the brief requires"},
            "notes": f"OBS quartile cut points on domain: {[float(np.trunc(c * 1e4) / 1e4) for c in cuts]} (truncated)"}

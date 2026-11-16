"""M1 builder, e1d9d038bb (UHR and DKD in type 2 diabetes, Medicine 2026, doi [DOI withheld]).
Spec: ln-UHR Q4 vs Q1, DKD (original) = UACR >= 30 and/or eGFR < 60 (CKD-EPI 2009 without race multiplier),
NHANES 2009-2018, Model 3 (primary), WTMEC2YR/5 with SDMVPSU and SDMVSTRA.
Sources in fulltext_txt/e1d9d038bb.txt: design and weights (2.1, line 13), exclusions (2.2, line 16), UHR, diabetes
and DKD (2.3, lines 19-21), covariates (2.4), models (2.5, line 29), Table 1 (variable coding), Table 2.
"""
import numpy as np
import pandas as pd

REC_ID = "e1d9d038bb"
CYCLES = ["2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def _drinks_per_day(df, lib):
    """Average drinks per drinking day in the past 12 months (ALQ130), 0 for people who did not drink in the past
    12 months (2009-2016: ALQ110 = 2 never 12 drinks, or ALQ120Q = 0; 2017-2018: ALQ111 = 2 never, or ALQ121 = 0)."""
    s = lib.stems("ALQ")
    a101, a110 = lib.raw(df, s, "ALQ101"), lib.raw(df, s, "ALQ110")
    a120 = lib.raw(df, s, "ALQ120Q")
    a111, a121 = lib.raw(df, s, "ALQ111"), lib.raw(df, s, "ALQ121")
    a130 = lib.raw(df, s, "ALQ130")
    dpd = a130.where(a130.between(0, 100))
    none = (((a101 == 2) & (a110 == 2)) | (a120 == 0) | (a111 == 2) | (a121 == 0))
    return pd.Series(np.where(none, 0.0, dpd), index=df.index)


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    uhr = lib.INDEX["UHR_x100"][0](df)
    df["ln_uhr"] = np.log(uhr.where(uhr > 0))
    # diabetes: prior diagnosis, FPG >= 7.0 mmol/L, HbA1c >= 6.5 %, diabetes medication or insulin
    diq070 = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    df["diab"] = ((df.DIQ010 == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5) | (diq070 == 1) | (df.DIQ050 == 1)).astype(float)
    # eGFR: CKD-EPI 2009 sex-specific, race multiplier not applied to anyone
    df["egfr_norace"] = lib.nc.egfr_ckdepi2009(df.scr, df.age, df.female == 1, pd.Series(False, index=df.index))
    e = df.egfr_norace
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())
    cancer = lib.raw(df, s("MCQ"), "MCQ220")
    df["smk100"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    df["drink"] = _drinks_per_day(df, lib)
    dom = ((df.age >= 20) & (df.pregnant != 1) & df.ln_uhr.notna() & df.uacr.notna() & df.egfr_norace.notna()
           & (cancer != 1) & (df.diab == 1))
    qnames, cuts = lib.quantile_dummies(df, "ln_uhr", dom, k=4, prefix="lnuhr_q")
    terms = qnames + ["age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["pir", "smk100", "drink", "bmi", "hba1c", "hypertension_any_140"]
    df["w"] = lib.pooled_weight(df, "mec")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "lnuhr_q4", "weight": "w",
            "readings": {
                "which model": "reading_1: Model 3, Q4 vs Q1 OR 2.43",
                "exposure scale": "reading_1: ln-UHR quartiles, Q4 vs Q1; unweighted quartiles in the domain",
                "T2DM": "the paper's diabetes definition (diagnosis, FPG >= 126 mg/dL, HbA1c >= 6.5, pills DIQ070, insulin DIQ050); no type 1 separation stated",
                "cancer": "MCQ220 = 1 excluded; unknown or missing MCQ220 kept",
                "eGFR equation": "CKD-EPI 2009 without the race multiplier, as written",
                "race": "5 groups of Table 1 (RIDRETH1)",
                "education": "3 groups of Table 1 (educ3)",
                "smoking": "smoked at least 100 cigarettes in life (SMQ020), the Table 1 coding",
                "drinking": "Table 1 reports drinking as a mean: average drinks per drinking day past 12 months (ALQ130), 0 for those not drinking in the past 12 months",
                "age, PIR, BMI, HbA1c": "continuous (Methods give no categories; Table 1 groups are descriptive)",
                "hypertension": "not defined in the paper: common convention diagnosis, current medication or mean BP >= 140/90 (hypertension_any_140)",
                "weight": "WTMEC2YR/5 (pooled MEC weight)",
                "missing data": "complete case (paper's primary analysis; multiple imputation only as a check)",
            },
            "notes": f"ln-UHR quartile cut points in the domain: {[int(c * 1000) / 1000 for c in cuts]}"}

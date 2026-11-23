"""M1 builder, 450ed6ebb6 (remnant cholesterol and DKD; source withheld).
Spec: RC quartile Q4 vs Q1, DKD (diabetes and (ACR >= 30 or eGFR <= 60)), Model 3 (11 covariates),
NHANES 2015-2016 and 2017-March 2020, weighted.
Sources in the paper (fulltext_txt/450ed6ebb6.txt): population line 12, RC line 15, quartiles line 18, DKD lines 21-24,
covariates line 27, weights and models lines 30-31.
"""
import numpy as np
import pandas as pd

REC_ID = "450ed6ebb6"
CYCLES = ["2015-2016", "2017-2020"]


def _miss(v, codes=(7, 9, 77, 99)):
    return v.where(~v.isin(codes))


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["RC"][0](df)
    # diabetes (ADA): self-reported diagnosis, current antidiabetic medication (insulin or pills), HbA1c >= 6.5 %
    pos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (df.DIQ070 == 1) | (df.hba1c >= 6.5)
    df["dm"] = pd.Series(np.where(pos, 1.0, np.where(df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna(), 0.0, np.nan)), index=df.index)
    e = df.egfr_ckdepi2009
    kid = (df.uacr >= 30) | (e <= 60)
    df["y"] = ((df.dm == 1) & kid).astype(float).where(df.dm.notna() & df.uacr.notna() & e.notna())

    # marital: married (incl. living with partner, the 2017-2020 code DMDMARTZ = 1 does not separate them) vs other
    m1 = _miss(lib.raw(df, s("DEMO"), "DMDMARTL"))
    m2 = _miss(lib.raw(df, s("DEMO"), "DMDMARTZ"))
    df["married"] = np.where(df.cycle == "2015-2016", m1.isin([1, 6]).astype(float).where(m1.isin([1, 2, 3, 4, 5, 6])),
                             (m2 == 1).astype(float).where(m2.isin([1, 2, 3])))
    # education <= high school diploma vs > high school
    df["educ_gt_hs"] = (df.educ3 == "gt_hs").astype(float).where(df.educ3.notna())
    # PIR < 1 vs >= 1
    df["pir_lt1"] = (df.pir < 1).astype(float).where(df.pir.notna())
    # smoking yes/no: smoked >= 100 cigarettes in life
    df["smoke_yes"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    # drinking yes/no: any alcohol in the past 12 months
    # 2015-2016: ALQ120Q > 0 (0 = none; never 12 drinks in a year nor lifetime = no); 2017-2020: ALQ121 1-10 (0 = none), ALQ111 = 2 no
    a101 = _miss(lib.raw(df, s("ALQ"), "ALQ101")); a110 = _miss(lib.raw(df, s("ALQ"), "ALQ110"))
    q120 = lib.raw(df, s("ALQ"), "ALQ120Q"); q120 = q120.where(~q120.isin([777, 999]))
    a111 = _miss(lib.raw(df, s("ALQ"), "ALQ111")); a121 = _miss(lib.raw(df, s("ALQ"), "ALQ121"))
    old = pd.Series(np.where((q120 > 0), 1.0, np.where((q120 == 0) | ((a101 == 2) & (a110 == 2)), 0.0, np.nan)), index=df.index)
    new = pd.Series(np.where(a121.between(1, 10), 1.0, np.where((a121 == 0) | (a111 == 2), 0.0, np.nan)), index=df.index)
    df["drink_yes"] = np.where(df.cycle == "2015-2016", old, new)
    # hypertension: self-reported diagnosis, or SBP >= 140 and/or DBP >= 90 (no medication criterion)
    pos = (df.BPQ020 == 1) | (df.sbp >= 140) | (df.dbp >= 90)
    df["htn"] = pd.Series(np.where(pos, 1.0, np.where(df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna()), 0.0, np.nan)), index=df.index)
    # diabetes duration (years) = age at survey - age at diagnosis (DID040; 666 = under 1 year old = 0);
    # 0 for participants without diagnosed diabetes
    did = lib.raw(df, s("DIQ"), "DID040")
    age_dx = did.where(~did.isin([777, 999])).replace(666, 0)
    dur = (df.age - age_dx).clip(lower=0)
    df["dm_dur"] = np.where(df.DIQ010 == 1, dur, np.where(df.DIQ010.isin([2, 3]) | (df.dm == 0), 0.0, np.nan))
    df["dm_dur"] = df["dm_dur"].astype(float)

    dom = ((df.age >= 20) & (df.pregnant != 1) & e.notna() & df.uacr.notna() & df.albumin_serum.notna() & df.scr.notna()
           & df.dm.notna() & df.crp_mgL.notna() & df.tc.notna() & df.hdl.notna() & df.ldl.notna() & df.x.notna())
    qn, cuts = lib.quantile_dummies(df, "x", dom, k=4)

    terms = qn + ["age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["educ_gt_hs", "married", "pir_lt1", "smoke_yes", "drink_yes", "bmi", "htn", "dm_dur"]

    df["w"] = lib.pooled_weight(df, "fasting")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x_q4", "weight": "w",
            "readings": {
                "main estimate": "reading_1 index RC, Q4 vs Q1, Model 3 OR 2.73",
                "DKD definition": "reading_1 Methods: diabetes and (ACR >= 30 or eGFR <= 60, 'less than or equal' as written)",
                "eGFR equation": "CKD-EPI 2009 with race (the formula printed in the paper has the 1.159 black term)",
                "population": "all adults >= 20, not pregnant (the paper's population is not restricted to diabetes: 1,014 DKD of 5,943); recon CFG diabetes restriction not used",
                "domain completeness": "eGFR, ACR, serum albumin (ALB read as serum albumin), Scr, diabetes status, hs-CRP, TC, HDL-C, LDL-C non-missing",
                "quartiles": "unweighted RC quartiles cut on the domain ('based on their distribution in the study population')",
                "survey weighting": "weights not named; RC needs LDL-C (fasting subsample): pooled fasting weight WTSAF2YR / WTSAFPRP",
                "diabetes": "DIQ010 = 1, insulin DIQ050 or pills DIQ070, HbA1c >= 6.5 (fasting glucose not in the definition)",
                "race": "RIDRETH1 5 groups as listed",
                "marital status": "married or living with partner vs other (2017-2020 DMDMARTZ cannot separate them)",
                "education": "<= high school diploma (DMDEDUC2 1-3) vs > high school (4-5)",
                "PIR": "binary < 1 vs >= 1",
                "smoking": "yes/no not defined: smoked >= 100 cigarettes (SMQ020 = 1)",
                "drinking": "yes/no not defined: any alcohol in the past 12 months",
                "hypertension": "as written: diagnosis or SBP >= 140 and/or DBP >= 90",
                "diabetes duration": "age minus DID040 for diagnosed diabetes; 0 for participants without diagnosed diabetes (undiagnosed HbA1c-only diabetes included)",
                "BMI, age": "continuous",
                "missing data": "paper imputed covariates (na.roughfix); M1 fit is complete case (automatic)",
            },
            "notes": "builder written before first execution"}

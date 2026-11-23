"""M1 builder, b966d3cfdc (OBS and CKD in diabetes; source withheld).
Spec: OBS quartile Q4 vs Q1, CKD (eGFR < 60 or ACR >= 30), Model 3 (13 covariates), NHANES 2007-2018, weighted.
Sources in the paper (fulltext_txt/b966d3cfdc.txt): population p. 7 (lines 573-581), OBS p. 7 (lines 581-714),
diabetes and CKD (lines 716-721), covariates (lines 726-734), statistics p. 8 (quartiles Q1 <= P25 ..., Model 3).
OBS from diet_indices_i139.DIET['b966d3cfdc'] reading 1 (Table 5 fixed cut points, day-1 recall, leisure PA).
"""
import numpy as np
import pandas as pd

REC_ID = "b966d3cfdc"
CYCLES = ["2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]
OGTT = {"2007-2008": "OGTT_E", "2009-2010": "OGTT_F", "2011-2012": "OGTT_G", "2013-2014": "OGTT_H", "2015-2016": "OGTT_I"}


def _miss(v, codes=(7, 9, 77, 99)):
    return v.where(~v.isin(codes))


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.DIET[REC_ID](df, reading=1)
    e = df.egfr_ckdepi2009
    df["y"] = ((e < 60) | (df.uacr >= 30)).astype(float).where(df.uacr.notna() & e.notna())

    # diabetes (paper's definition): diagnosis, FPG >= 7.0 mmol/L, HbA1c >= 6.5 %, OGTT 2 h >= 11.1 mmol/L, insulin or pills
    ogtt = lib.raw(df, OGTT, "LBXGLT")
    df["dm"] = ((df.DIQ010 == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5) | (ogtt >= 200) | (df.DIQ050 == 1)
                | (df.DIQ070 == 1)).astype(float)
    # hypertension: SBP >= 140, DBP >= 90, or antihypertensive drugs (no self-reported diagnosis in the definition)
    med = df.BPQ_med_now == 1
    bp_hi = (df.sbp >= 140) | (df.dbp >= 90)
    df["htn"] = pd.Series(np.where(med | bp_hi, 1.0, np.where(df.sbp.notna() & df.dbp.notna(), 0.0, np.nan)), index=df.index)
    # smoking yes/no: smoked >= 100 cigarettes in life (SMQ020)
    df["smoke_yes"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    # BMI <= 30 vs > 30
    df["bmi_gt30"] = (df.bmi > 30).astype(float).where(df.bmi.notna())

    dom = (df.dm == 1) & (df.age >= 18) & df.x.notna() & df.uacr.notna() & e.notna()
    qn, cuts = lib.quantile_dummies(df, "x", dom, k=4)

    terms = qn + ["age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["htn", "smoke_yes", "bmi_gt30", "u_cr", "uric_acid", "albumin_serum", "bilirubin_umol", "tc", "bun"]

    df["w"] = lib.pooled_weight(df, "mec")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x_q4", "weight": "w",
            "readings": {
                "main estimate": "reading_1 Model 3 Q4 vs Q1 OR 0.69",
                "cycles": "reading_1 six two-year cycles 2007-2018",
                "CKD criteria": "reading_1 either criterion: eGFR < 60 or ACR >= 30",
                "eGFR equation": "CKD-EPI version not stated: 2009 with race",
                "OBS": "DIET['b966d3cfdc'] reading 1 (Table 5 sex-specific cut points, day-1 recall, leisure PA, all 17 components)",
                "quartiles": "unweighted quartiles of OBS cut on the domain; Q1 <= P25 (value at a cut point goes to the lower group, as written)",
                "survey weighting": "'recommended sample weight calculation', variable not named: pooled MEC weight (smallest lib kind; day-1 dietary weight is not a lib kind)",
                "diabetes": "as written; OGTT only 2007-2016 (no OGTT file 2017-2018); pills DIQ070 absent 2007-2008",
                "pregnancy": "not excluded (the paper lists no such exclusion)",
                "race": "5 groups as listed (RIDRETH1, Other incl. multiracial)",
                "education": "3 levels as listed (DMDEDUC2 1-2 / 3 / 4-5)",
                "hypertension": "as written: SBP >= 140, DBP >= 90 or antihypertensive drugs (BPQ050A); self-report of diagnosis not used",
                "smoking status": "yes/no not defined; read as smoked >= 100 cigarettes in life (SMQ020 = 1)",
                "BMI": "binary <= 30 vs > 30 as written",
                "laboratory covariates": "UCR urine creatinine URXUCR, UA, serum ALB, TBIL LBDSTBSI, TC, BUN entered linearly",
            },
            "notes": "builder written before first execution"}

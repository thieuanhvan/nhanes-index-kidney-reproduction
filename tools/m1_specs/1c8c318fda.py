"""M1 builder, 1c8c318fda (Ren Fail 2025, doi [DOI withheld]). UHR per SD, CKD, Model 3, NHANES part.
Spec: NHANES 1999-2010, age >= 45, CKD = eGFR < 60 (2021 CKD-EPI), weighted ('proper sample weights', not named),
Table 3 NHANES Model 3 OR 1.81 per SD.
Sources (fulltext_txt/1c8c318fda.txt): exclusions line 18, outcome line 26, covariates lines 29-30, models line 34,
Table 1 (NHANES) lines 99-156, Table 3 footnote line 221.
"""
import numpy as np, pandas as pd

REC_ID = "1c8c318fda"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010"]


def build(lib):
    df = lib.table(CYCLES)
    df["uhr"] = lib.INDEX["UHR"][0](df)
    e = df.egfr_ckdepi2021
    df["y"] = (e < 60).astype(float).where(e.notna())
    dom = (df.age >= 45) & df.uhr.notna() & df.scr.notna()
    df["x"] = lib.sd_scale(df, "uhr", dom)
    # marital status: married vs others (unmarried, cohabiting, divorced, separated, widowed)
    mar = lib.raw(df, lib.stems("DEMO"), "DMDMARTL")
    df["married"] = (mar == 1).astype(float).where(mar.isin([1, 2, 3, 4, 5, 6]))
    # smoking: ever (current or former) vs never
    df["ever_smoke"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    # drinking: ever (current or former) vs never; never = fewer than 12 drinks in any year and in lifetime
    s = lib.stems("ALQ")
    yr12 = lib.raw(df, s, ["ALQ101", "ALQ100", "ALD100"])
    life12 = lib.raw(df, s, "ALQ110")
    ever = (yr12 == 1) | (life12 == 1)
    never = (life12 == 2) & ((yr12 == 2) | yr12.isna())
    df["ever_drink"] = pd.Series(np.where(ever, 1.0, np.where(never, 0.0, np.nan)), index=df.index)
    # lipid-lowering drugs: taking prescribed medicine for high cholesterol (BPQ100D); not asked (skip) = 0
    l = lib.raw(df, lib.stems("BPQ"), "BPQ100D")
    rec = df.BPQ020.isin([1, 2])
    df["lipid_drug"] = pd.Series(np.where(l == 1, 1.0, np.where(rec & ~l.isin([7, 9]), 0.0, np.nan)), index=df.index)
    terms = ["x", "female", "age"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["married", "bmi", "tg_serum", "tc", "cvd_history", "diabetes_any", "hypertension_any_140", "crp_mgL",
              "lipid_drug", "ever_smoke", "ever_drink"]
    df["w"] = lib.pooled_weight(df, "mec")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "cohort (spec ambiguity)": "reading_1: NHANES part, Model 3 per SD OR 1.81, N 10,968",
                "Model 3 covariates (spec ambiguity)": "reading_1: Model 2 covariates (age, sex, race, education, marital status, BMI) plus the nine added",
                "age": "'age < 45' excluded read as age >= 45 kept",
                "per SD": "UHR (SUA mg/dL / HDL-C mg/dL) / SD over the domain rows (unweighted); the x100 'UHR%' label does not change the per-SD OR",
                "eGFR": "CKD-EPI 2021 race-free (stated), calibrated creatinine",
                "weight": "'proper sample weights' not named; MEC exam weight (all inputs are MEC labs), WTMEC4YR for 1999-2002",
                "triglycerides": "serum biochemistry TG (BIOPRO LBXSTR), not the fasting-subsample TG: the paper excludes only missing HDL-C, SUA, Scr and key covariates (fasting status is not an exclusion) and the method of TG is not stated",
                "race": "RIDRETH1 5 groups as in Table 1",
                "education": "less than high school / high school or equivalent / above high school (educ3)",
                "marital status": "DMDMARTL 1 married vs 2-6 others (living with partner in others, as stated); 77/99 missing",
                "BMI, TC, TG, CRP, age": "continuous (standardisation of covariates does not change the exposure OR)",
                "CVD": "self-report of CHF, CHD, angina, MI or stroke (cvd_history); the 'medication history' part is not specified",
                "DM": "fasting glucose >= 126 or HbA1c >= 6.5 or DIQ010 = 1 or insulin/pills (diabetes_any)",
                "hypertension": "mean SBP >= 140 or DBP >= 90, BPQ020 = 1 or BPQ050A = 1 (hypertension_any_140)",
                "CRP": "CRP mg/L (LBXCRP x 10)",
                "lipid-lowering drugs": "BPQ100D = 1 'taking prescriptions for hypercholesterolemia'; skipped = 0",
                "smoking": "ever (SMQ020 = 1) vs never",
                "drinking": "ever (12+ drinks in any year or in lifetime) vs never",
            },
            "notes": "written before first execution; no post-run changes"}

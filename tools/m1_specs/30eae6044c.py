"""M1 builder, 30eae6044c (ABSI and DKD in diabetes, Ren Fail 2025, doi [DOI withheld]).
Spec: ABSI x100 per 1 unit, DKD (diabetes plus eGFR2021 < 60 or ACR >= 30), Model 3, NHANES 2005-2006 to 2015-2016, unweighted.
Sources in fulltext_txt/30eae6044c.txt: exclusions line 12, ABSI line 15, outcome and diabetes line 18, covariates line 21,
models and mean imputation line 24, Table 2 footnote line 132.
"""
import numpy as np
import pandas as pd

REC_ID = "30eae6044c"
CYCLES = ["2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016"]


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["ABSI_x100"][0](df)
    df["y"] = ((df.egfr_ckdepi2021 < 60) | (df.uacr >= 30)).astype(float).where(df.uacr.notna() & df.egfr_ckdepi2021.notna())
    s = lib.stems
    # diabetes: diagnosis, insulin or pills, OGTT 2 h >= 200 mg/dL, FPG >= 126, HbA1c >= 6.5
    pills = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    ogtt = lib.raw(df, s("OGTT"), "LBXGLT")
    diab = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (pills == 1) | (ogtt >= 200) | (df.glucose >= 126) | (df.hba1c >= 6.5)
    cancer = lib.raw(df, s("MCQ"), "MCQ220") == 1
    # smoking: never/former vs current (SMQ040 1-2)
    df["smk_current"] = (df.smoking == "current").astype(float).where(df.smoking.notna())
    # drinking: never (< 12 drinks lifetime, ALQ110 == 2) vs drinker (ALQ101 == 1 or ALQ110 == 1)
    a101 = lib.raw(df, s("ALQ"), "ALQ101")
    a110 = lib.raw(df, s("ALQ"), "ALQ110")
    drk = pd.Series(np.nan, index=df.index)
    drk[(a101 == 1) | (a110 == 1)] = 1.0
    drk[(a110 == 2) & (a101 != 1)] = 0.0
    df["drinker"] = drk
    # hypertension: prior diagnosis or current antihypertensive medication
    htn = (df.BPQ020 == 1) | (df.BPQ_med_now == 1)
    df["htn"] = htn.astype(float).where(htn | df.BPQ020.isin([1, 2]))
    dom = ((df.age >= 20) & (df.age <= 80) & ~cancer & (df.pregnant != 1) & diab & df.y.notna() & df.x.notna()
           & df.smk_current.notna() & df.drinker.notna() & df.educ3.notna())
    # the paper's mean imputation is not used (prereg section 4: missing stays missing, complete case)
    df["pir3"] = pd.Series(np.select([df.pir <= 1.3, df.pir <= 3.5, df.pir > 3.5], ["low", "mid", "high"], default=""),
                           index=df.index, dtype=object).replace("", np.nan)
    terms = ["x", "female", "age"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "pir3", "low")
    terms += ["smk_current", "drinker", "htn", "cvd_history", "alt", "ast", "tc", "hdl"]
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": None,
            "readings": {
                "survey weighting": "none (reading_1)",
                "index components": "formula inputs waist, BMI, height (reading_1), INDEX ABSI_x100",
                "age": "20 to 80 inclusive (80 is the 80+ top code)",
                "diabetes": "DIQ010 == 1, insulin DIQ050, pills DIQ070/DID070, OGTT LBXGLT >= 200, FPG >= 126, HbA1c >= 6.5",
                "cancer": "MCQ220 == 1 excluded",
                "eGFR": "CKD-EPI 2021 as stated; CKD = eGFR < 60 or ACR >= 30",
                "smoking": "binary current (SMQ040 every day or some days) vs never/former",
                "drinking": "never (ALQ110 == 2, < 12 drinks lifetime) vs drinker (ALQ101 == 1 or ALQ110 == 1)",
                "hypertension": "BPQ020 == 1 or BPQ050A == 1 (no measured BP, as written)",
                "CVD": "cvd_history MCQ160B-F",
                "race, education": "RIDRETH1 5 groups; DMDEDUC2 3 groups",
                "PIR": "3 categories <= 1.3, 1.3-3.5, > 3.5",
                "mean imputation": "not used: prereg section 4 complete case overrides the paper's mean imputation",
                "BMI, WC, HbA1c": "listed as covariates in Methods but not in the Model 3 list; not used"},
            "notes": "builder by agent 15/11. Post-run change required by prereg section 4 missing-data rule, not by n: removed the mean imputation of ALT, AST, TC, HDL, PIR; missing stays NaN and complete case drops it"}

"""M1 builder, ea03123ab4 (METS-VF and albuminuria, Front Endocrinol 2025, doi [DOI withheld]).
Spec: METS-VF continuous per 1 unit, Model 3, albuminuria UACR >= 30, NHANES 1999-2018.
Sources (fulltext_txt/ea03123ab4.txt): population and formula (Methods), covariates section, statistical analysis
(weights stated for descriptive statistics only; Model 3 list), Table 1 coding, Table 3.
"""
import numpy as np
import pandas as pd

REC_ID = "ea03123ab4"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010",
          "2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def labels(index, pairs):
    out = pd.Series(np.nan, index=index, dtype=object)
    for m, lab in reversed(pairs):
        out[m.fillna(False).to_numpy(bool)] = lab
    return out


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["METS-VF_3.239"][0](df)
    df["y"] = (df.uacr >= 30).astype(float).where(df.uacr.notna())
    dom = (df.age >= 20) & df.x.notna() & df.uacr.notna()
    df["age60"] = flag(df.age >= 60, df.age.notna())
    mar = lib.raw(df, lib.stems("DEMO"), "DMDMARTL")
    df["married"] = flag(mar == 1, mar.isin([1, 2, 3, 4, 5, 6]))
    smq = lib.raw(df, lib.stems("SMQ"), "SMQ020")
    df["smoke_hist"] = flag(smq == 1, smq.isin([1, 2]))
    df["htn"] = df["hypertension_any_140"]
    pills = lib.raw(df, lib.stems("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, lib.stems("DIQ"), "DIQ050")
    df["dm"] = flag((df.glucose >= 126) | (df.hba1c >= 6.5) | (df.dx_diabetes == 1) | (pills == 1) | (ins == 1),
                    df.dx_diabetes.notna() | df.hba1c.notna() | df.glucose.notna())
    df["bmi3"] = labels(df.index, [(df.bmi < 25, "lt25"), (df.bmi < 30, "25to30"), (df.bmi >= 30, "ge30")])
    terms = ["x", "age60", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["married", "pir"]
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["smoke_hist", "htn", "dm", "cvd_history"]
    terms += lib.dummies(df, "bmi3", "lt25")
    terms += ["tc", "ldl", "egfr_ckdepi2009"]
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": None,
            "readings": {
                "survey weights in regression (spec ambiguity)": "reading_1: unweighted regression (weights stated for descriptive statistics only)",
                "index components (spec ambiguity)": "reading_1: base variables FPG, TG, BMI, HDL-C, WC, height, sex, age (recon key METS-VF_3.239)",
                "UACR threshold": ">= 30 mg/g",
                "age": "age group < 60 / >= 60 (Model 3 list)",
                "ethnicity": "RIDRETH1 5 groups, reference NH white",
                "marital status": "Table 1 'Married' yes/no: DMDMARTL = 1 vs 2-6 (living with partner counted as not married)",
                "PIR": "continuous (Table 1 median, no categories stated)",
                "education": "3 groups: less than high school, high school, some college or above",
                "smoking history": "current or former smoker = SMQ020 = 1 (at least 100 cigarettes) vs never",
                "hypertension": "SBP >= 140, DBP >= 90 (mean of readings), BPQ020 = 1, or BPQ050A = 1",
                "diabetes": "fasting glucose >= 126, HbA1c >= 6.5, DIQ010 = 1, pills (DIQ070; DID070 in 2005-2008) or insulin (DIQ050)",
                "CVDs": "cvd_history (MCQ160B-F)",
                "BMI group": "< 25 (reference), 25 to < 30, >= 30",
                "TC, LDL-c": "continuous mg/dL",
                "eGFR": "CKD-EPI 2009 with race (version not stated), continuous",
                "missing covariates": "paper imputed with random forest; here complete case on model terms"},
            "notes": "domain needs the fasting subsample (FPG and TG in METS-IR); regression unweighted per reading_1"}

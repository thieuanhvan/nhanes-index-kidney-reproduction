"""M1 builder, c62af9091f (CALLY index and CKD, US NHANES 1999-2010 arm; doi [DOI withheld]).
Spec: CALLY continuous per 1 unit, Model 3, CKD = eGFR < 60 and/or UACR >= 30.
Sources (fulltext_txt/c62af9091f.txt): exclusions line 15, outcome line 18, index line 20, covariates line 25,
models line 29, Table 1 lines 106-138, Table 3 lines 185-194. No survey weighting is mentioned anywhere.
"""
import numpy as np
import pandas as pd

REC_ID = "c62af9091f"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010"]


def flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def labels(index, pairs):
    out = pd.Series(np.nan, index=index, dtype=object)
    for m, lab in reversed(pairs):
        out[m.fillna(False).to_numpy(bool)] = lab
    return out


def build(lib):
    df = lib.table(CYCLES)
    # CALLY = lymphocytes [10^9/L] x albumin [g/dL] / (CRP [mg/dL] x 10); crp_mgL / 10 = mg/dL
    df["x"] = lib.INDEX["CALLY_gdL"][0](df)
    e = df.egfr_ckdepi2009
    df["y"] = ((e < 60) | (df.uacr >= 30)).astype(float).where(e.notna() & df.uacr.notna())
    dom = (df.age >= 20) & (df.pregnant != 1) & df.scr.notna() & df.uacr.notna() & df.x.notna()
    # race: 4 groups (Mexican American, NH white, NH black, other incl. other Hispanic)
    df["race4"] = df.race1.map({"mexican_american": "mex", "nh_white": "white", "nh_black": "black",
                                "other_hispanic": "other", "other_multi": "other"})
    mar = lib.raw(df, lib.stems("DEMO"), "DMDMARTL")
    df["cohab"] = flag(mar.isin([1, 6]), mar.isin([1, 2, 3, 4, 5, 6]))
    smq = lib.raw(df, lib.stems("SMQ"), "SMQ020")
    df["smoke_ever"] = flag(smq == 1, smq.isin([1, 2]))
    alq = lib.raw(df, lib.stems("ALQ"), ["ALQ101", "ALQ100", "ALD100"])
    df["alcohol"] = flag(alq == 1, alq.isin([1, 2]))
    df["htn"] = df["hypertension_any_140"]
    pills = lib.raw(df, lib.stems("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, lib.stems("DIQ"), "DIQ050")
    df["dm"] = flag((df.dx_diabetes == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5) | (pills == 1) | (ins == 1),
                    df.dx_diabetes.notna() | df.hba1c.notna() | df.glucose.notna())
    low_hdl = np.where(df.female == 1, df.hdl < 50, df.hdl < 40)
    pos = (df.tc >= 200) | (df.tg >= 150) | pd.Series(low_hdl, index=df.index) & df.hdl.notna() | (df.ldl >= 130)
    df["dyslip"] = flag(pos, df.tc.notna() & df.hdl.notna())
    df["obese28"] = flag(df.bmi >= 28, df.bmi.notna())
    terms = ["x", "female", "age"]
    terms += lib.dummies(df, "race4", "white")
    terms += ["cohab", "smoke_ever", "alcohol", "htn", "dm", "dyslip", "obese28"]
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": None,
            "readings": {
                "eGFR equation (spec ambiguity)": "reading_1: 'most recent iteration', version not stated, so CKD-EPI 2009 with race (brief rule 5)",
                "UACR threshold": ">= 30 mg/g",
                "survey weighting": "unweighted (no weighting mentioned; spec weight_variable null)",
                "age": "continuous (Table 1 median, no categories)",
                "race": "4 groups as listed: Mexican American, NH white (reference), NH black, other (other Hispanic and other/multiracial)",
                "marital status": "Table 1 'Cohabitation' vs 'Solitude': married or living with partner (DMDMARTL 1, 6) vs widowed, divorced, separated, never married",
                "smoking": "yes = smoked at least 100 cigarettes in life (SMQ020 = 1), Table 1 yes/no",
                "alcohol": "yes = at least 12 drinks in any one year (ALQ100 1999-2000, ALD100 2001-2002, ALQ101 2003-2010), Table 1 yes/no; definition not given in the paper",
                "hypertension": "BPQ020 = 1, SBP >= 140 or DBP >= 90 (mean of readings), or BPQ050A = 1",
                "diabetes": "DIQ010 = 1, fasting glucose >= 126, HbA1c >= 6.5, or pills (DIQ070; DID070 in 2005-2008) or insulin (DIQ050)",
                "dyslipidemia": "any of TC >= 200, TG >= 150, HDL < 50 (women) / < 40 (men), LDL >= 130; known when TC and HDL present (TG, LDL only in the fasting subsample; missing counts as not met)",
                "obesity": "BMI >= 28 (as stated)",
                "exclusion 4 (incomplete demographics or chronic diseases)": "left to complete case on model terms"},
            "notes": "CRP 1999-2010 from LAB11/L11_B/L11_C/CRP_D-F (LBXCRP mg/dL x 10 = mg/L in crp_mgL)"}

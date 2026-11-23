"""M1 builder, 6b7f67827b (NTI and CKD in hypertensive adults; source withheld).
Spec: NHANES part (reading_1), NTI per 1 unit, CKD = eGFR < 60 or UACR >= 30, Model 2, NHANES 1999-2000 to 2017-March 2020,
pooled fasting-subsample weights (2017-2020 as 3.2 years).
Sources in fulltext_txt/6b7f67827b.txt: population line 17, exposure line 20, outcome line 23, covariates lines 26-27,
models line 32, weights line 36, Table 2 (NHANES categories incl. "Unknown" levels) lines 135-181, Table 3 footnote line 203.
"""
import numpy as np
import pandas as pd

REC_ID = "6b7f67827b"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014",
          "2015-2016", "2017-2020"]
ALQ121_DAYS = {0: 0, 1: 365, 2: 286, 3: 182, 4: 104, 5: 52, 6: 30, 7: 12, 8: 9, 9: 4.5, 10: 1.5}


def _miss(v, codes):
    return v.where(~v.isin(codes))


def _drink_days(lib, df):
    """Drinking days in the past 12 months (ALQ120Q/U 1999-2016, ALQ121 2017+); 0 when skipped as non-drinker."""
    s = lib.stems("ALQ")
    q = _miss(lib.raw(df, s, "ALQ120Q"), [777, 999])
    u = _miss(lib.raw(df, s, "ALQ120U"), [7, 9])
    days = q * u.map({1: 52, 2: 12, 3: 1})
    days = days.where(q != 0, 0.0)
    q121 = _miss(lib.raw(df, s, "ALQ121"), [77, 99])
    days = days.fillna(q121.map(ALQ121_DAYS))
    skip = (lib.raw(df, s, "ALQ110") == 2) | (lib.raw(df, s, "ALQ111") == 2)
    return days.where(~(skip & days.isna()), 0.0)


def _cat(cond_list, labels, df):
    # post-run change (prereg section 4): no missing level, unmatched rows stay NaN
    out = pd.Series(np.select(cond_list, labels, default=""), index=df.index, dtype=object)
    return out.replace("", np.nan)


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["NTI"][0](df)
    df["y"] = ((df.egfr_ckdepi2009 < 60) | (df.uacr >= 30)).astype(float).where(df.uacr.notna() & df.egfr_ckdepi2009.notna())
    s = lib.stems
    # hypertension (population): SBP >= 140, DBP >= 90, self-reported diagnosis, or current antihypertensive medication
    htn = df.hypertension_any_140 == 1
    # dialysis (KIQ025, past 12 months; file absent in 1999-2000)
    dial = lib.raw(df, s("KIQ_U"), "KIQ025") == 1
    pos_comp = (df.neut_n > 0) & (df.lymph_n > 0) & (df.tg > 0) & (df.glucose > 0)
    # covariates as coded in the paper; the paper's Unknown levels are not used (prereg section 4)
    df["educ2"] = _cat([df.educ_raw.isin([1, 2]), df.educ_raw.isin([3, 4, 5])], ["below_hs", "hs_above"], df)
    mar = lib.raw(df, s("DEMO"), "DMDMARTL")
    marz = lib.raw(df, s("DEMO"), "DMDMARTZ")
    df["married"] = ((mar == 1) | (marz == 1)).astype(float).where(mar.isin([1, 2, 3, 4, 5, 6]) | marz.isin([1, 2, 3]))
    df["nhw"] = (df.race1 == "nh_white").astype(float).where(df.race1.notna())
    df["smk"] = df.smoking.astype(object)
    dd = _drink_days(lib, df)
    df["drink_cur"] = (dd > 0).astype(float).where(dd.notna())
    pills = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    dpos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (pills == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5)
    dobs = df.DIQ010.isin([1, 2, 3]) | df.glucose.notna() | df.hba1c.notna()
    df["diab"] = _cat([dpos, dobs], ["yes", "no"], df)
    hf = lib.raw(df, s("MCQ"), "MCQ160B")
    df["hf"] = _cat([hf == 1, hf == 2], ["yes", "no"], df)
    terms = ["x", "age", "female", "nhw"]
    terms += lib.dummies(df, "educ2", "hs_above")
    terms += ["married"]
    terms += lib.dummies(df, "smk", "never")
    terms += ["drink_cur", "bmi", "sbp"]
    terms += lib.dummies(df, "diab", "no")
    terms += lib.dummies(df, "hf", "no")
    df["w"] = lib.pooled_weight(df, "fasting")
    dom = (df.age >= 20) & htn & pos_comp & df.x.notna() & ~dial
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "cohort": "NHANES part (reading_1)",
                "eGFR equation": "not stated for NHANES (reading_1); CKD-EPI 2009 with race per brief rule 5",
                "UACR threshold": ">= 30 mg/g as written",
                "hypertension": "hypertension_any_140: BPQ020, BPQ050A current medication, mean SBP >= 140 or DBP >= 90",
                "valid positive NTI components": "neutrophils, lymphocytes, TG, fasting glucose all > 0",
                "dialysis": "KIQ025 == 1 excluded (question absent in 1999-2000)",
                "education": "below high school (DMDEDUC2 1-2) / high school or above (3-5); missing stays missing (prereg section 4)",
                "marital": "married (DMDMARTL 1; 2017-2020 DMDMARTZ 1 includes living with partner) vs other; missing stays missing (prereg section 4)",
                "race": "non-Hispanic White vs other (RIDRETH1)",
                "smoking": "never/former/current from SMQ020/SMQ040; missing stays missing (prereg section 4)",
                "drinking": "current alcohol use = any drinking day in past 12 months (ALQ120Q > 0; ALQ121 1-10 in 2017-2020); skipped non-drinkers non-current; missing stays missing (prereg section 4)",
                "diabetes": "FPG >= 126, HbA1c >= 6.5, DIQ010 == 1, insulin DIQ050, pills DIQ070/DID070; missing when nothing observed",
                "heart failure": "MCQ160B yes/no; missing stays missing",
                "BMI, SBP": "continuous, complete case (paper excluded residual missing BMI or SBP)",
                "weight": "pooled fasting weight, 4-year weights for 1999-2002 per NCHS rule, 2017-2020 as 3.2 years"},
            "notes": "builder by agent 15/11; first version kept the paper's Unknown levels (Table 2). "
                     "Post-run technical fix (first m1_check): term diab__unknown was constant (all 0) in the domain, dropped; no reading changed (now moot: no unknown level exists). "
                     "Post-run change required by prereg section 4 missing-data rule, not by n: removed the Unknown levels (education, smoking, diabetes, heart failure) and the missing-as-other (marital) and missing-as-non-current (drinking) rules; missing stays NaN and complete case drops it"}

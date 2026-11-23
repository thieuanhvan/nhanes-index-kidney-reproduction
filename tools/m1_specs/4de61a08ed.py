"""M1 builder, 4de61a08ed (HRR and albuminuria; source withheld).
Spec: two-segment logistic model with the paper's inflection point 1.213; target = per 0.1 unit of HRR below 1.213,
albuminuria = UACR > 30 mg/g, NHANES 2007-2018, Table 3 adjustment list, unweighted (no weighting reported).
Sources in fulltext_txt/4de61a08ed.txt: population (line 12), HRR and outcome (lines 15-16), covariate definitions
(lines 19-23), model (Statistical analysis line 26; Results line 38), Table 3 footnote.
"""
import numpy as np
import pandas as pd

REC_ID = "4de61a08ed"
CYCLES = ["2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]
K = 1.213
# ALQ121 (2017-2018) frequency codes -> drinking days per year (code midpoints)
ALQ121_DAYS = {0: 0.0, 1: 365.0, 2: 286.0, 3: 182.0, 4: 104.0, 5: 52.0, 6: 30.0, 7: 12.0, 8: 9.0, 9: 4.5, 10: 1.5}


def _drinks_past_year(df, lib):
    """Number of alcoholic drinks in the past 12 months = drinking days x drinks per drinking day (ALQ130);
    0 for people with no drinking in the past 12 months or never 12 drinks."""
    s = lib.stems("ALQ")
    a101, a110 = lib.raw(df, s, "ALQ101"), lib.raw(df, s, "ALQ110")
    q, u = lib.raw(df, s, "ALQ120Q"), lib.raw(df, s, "ALQ120U")
    a111, a121 = lib.raw(df, s, "ALQ111"), lib.raw(df, s, "ALQ121")
    a130 = lib.raw(df, s, "ALQ130")
    dpd = a130.where(a130.between(0, 100))
    days_old = q.where(q.between(0, 366)) * u.map({1: 52.0, 2: 12.0, 3: 1.0})
    days_old = days_old.where(q != 0, 0.0)
    days_new = a121.map(ALQ121_DAYS)
    days = days_old.where(days_old.notna(), days_new)
    tot = days * dpd
    tot = tot.where(days != 0, 0.0)
    none = ((a101 == 2) & (a110 == 2)) | (a111 == 2)
    return pd.Series(np.where(none, 0.0, tot), index=df.index)


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["hrr"] = lib.INDEX["HRR"][0](df).round(2)   # "rounded to two decimal places"
    df["seg_lo"] = np.minimum(df.hrr, K) * 10          # per 0.1 unit below the inflection point
    df["seg_hi"] = np.maximum(df.hrr - K, 0) * 10      # per 0.1 unit above it
    df.loc[df.hrr.isna(), ["seg_lo", "seg_hi"]] = np.nan
    df["y"] = (df.uacr > 30).astype(float).where(df.uacr.notna())
    # smokers: at least 100 cigarettes in lifetime
    df["smoker"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    # drinkers: at least 12 alcoholic drinks in the past 12 months
    tot = _drinks_past_year(df, lib)
    df["drinker"] = (tot >= 12).astype(float).where(tot.notna())
    # hypertension: mean SBP >= 140 and/or DBP >= 90, or told 2 or more times (BPQ030), or told to take prescribed
    # antihypertensive medication (BPQ040A)
    b = s("BPQ")
    b030, b040 = lib.raw(df, b, "BPQ030"), lib.raw(df, b, "BPQ040A")
    known = (df.sbp.notna() & df.dbp.notna()) | df.BPQ020.isin([1, 2])
    df["htn"] = ((df.sbp >= 140) | (df.dbp >= 90) | (b030 == 1) | (b040 == 1)).astype(float).where(known)
    # diabetes: physician diagnosis, anti-diabetic pills, insulin, FBG >= 7.0 mmol/L or HbA1c >= 6.5 %
    diq070 = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    pos = (df.DIQ010 == 1) | (diq070 == 1) | (df.DIQ050 == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5)
    known = df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna()
    df["diab"] = pos.astype(float).where(known)
    terms = ["seg_lo", "seg_hi", "female", "age", "bmi", "pir", "tc", "tg", "alt", "ast", "bun", "uric_acid", "scr",
             "egfr_ckdepi2009"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["smoker", "drinker", "htn", "diab", "cvd_history"]
    terms += lib.dummies(df, "educ3", "lt_hs")
    dom = (df.age >= 20) & df.uacr.notna() & df.hrr.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "seg_lo", "weight": None,
            "readings": {
                "which segment": "reading_1: HRR < 1.213, OR 0.93 per 0.1 unit",
                "covariate set": "reading_1: Table 3 footnote (19 variables, no LDL-C, HDL-C, FBG, HbA1c)",
                "piecewise model": "two linear segments with the paper's inflection point 1.213 fixed: min(HRR, 1.213) and max(HRR - 1.213, 0), each x 10; target = lower segment",
                "HRR": "hemoglobin / RDW rounded to 2 decimals, as written",
                "outcome": "UACR > 30 mg/g as written",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated)",
                "TG": "triglycerides of the fasting lipid file (LBXTR, same file as LDL-C and FBG named in the Methods)",
                "race": "5 groups as written (RIDRETH1)",
                "education": "3 groups as written (educ3)",
                "smoking": "smoker = at least 100 cigarettes in lifetime, as written",
                "drinking": "drinker = at least 12 drinks in the past 12 months, as written: drinking days (ALQ120Q/U 2007-2016, ALQ121 code midpoints 2017-2018) x drinks per day (ALQ130) >= 12",
                "hypertension": "mean SBP >= 140 and/or DBP >= 90, or told 2 or more times (BPQ030), or told to take prescribed medicine (BPQ040A)",
                "diabetes": "diagnosis, pills DIQ070 (DID070 2007-2008), insulin DIQ050, FBG >= 126 mg/dL, HbA1c >= 6.5",
                "CVD": "CHD, angina, heart failure, MI or stroke (MCQ160B-F)",
                "continuous covariates": "age, BMI, PIR, TC, TG, ALT, AST, BUN, UA, Cr, eGFR continuous",
                "weight": "none (no survey weighting reported; spec weight_variable null)",
            },
            "notes": "segment coding gives the slope below and the slope above the knot as separate coefficients"}

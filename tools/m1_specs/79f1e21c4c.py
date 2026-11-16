"""M1 builder, 79f1e21c4c (Urol Nephrol 2024, doi [DOI withheld]). TyG per 1 unit, CKD, Model 3.
Spec: NHANES 1999-2020 (reading_1: 1999-2000 to 2015-2016 plus the 2017-March 2020 pre-pandemic file), CKD = eGFR < 60
(CKD-EPI 2009) or ACR >= 30 mg/g, weighted (variable not named), Table 2 Model 3 OR 4.03.
Sources (fulltext_txt/79f1e21c4c.txt): exclusions line 12, outcome line 16, covariates lines 19-20, models line 23,
Table 1 coding lines 102-156, Table 2 footnote line 179.
"""
import numpy as np, pandas as pd

REC_ID = "79f1e21c4c"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010",
          "2011-2012", "2013-2014", "2015-2016", "2017-2020"]


def _clean(s, bad):
    return s.where(~s.isin(bad))


def _alcohol_days(lib, df):
    """Days of alcohol consumption in the past year. 1999-2016: ALQ120Q x (52 week, 12 month, 1 year by ALQ120U);
    0 when ALQ120Q = 0 or when the person was skipped as a lifetime non-drinker (12-drink question(s) answered no).
    2017-2020: ALQ121 frequency category converted to days per year; 0 when ALQ111 = 2 (never had a drink)."""
    s = lib.stems("ALQ")
    q = _clean(lib.raw(df, s, "ALQ120Q"), [777, 999])
    u = _clean(lib.raw(df, s, "ALQ120U"), [7, 9])
    yr12 = _clean(lib.raw(df, s, ["ALQ101", "ALQ100", "ALD100"]), [7, 9])
    life12 = _clean(lib.raw(df, s, "ALQ110"), [7, 9])
    days = q * u.map({1: 52.0, 2: 12.0, 3: 1.0})
    days = days.where(~(q == 0), 0.0)
    never = (yr12 == 2) & (life12 == 2)
    days = days.where(~(never & q.isna()), 0.0)
    a121 = _clean(lib.raw(df, s, "ALQ121"), [77, 99])
    a111 = lib.raw(df, s, "ALQ111")
    per_year = {0: 0.0, 1: 365.0, 2: 312.0, 3: 182.0, 4: 104.0, 5: 52.0, 6: 30.0, 7: 12.0, 8: 9.0, 9: 4.5, 10: 1.5}
    p = df.cycle == "2017-2020"
    dp = a121.map(per_year)
    dp = dp.where(~((a111 == 2) & a121.isna()), 0.0)
    days[p] = dp[p]
    return days


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["TyG"][0](df)
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())
    # BMI as coded in Table 1: normal / overweight / obese
    df["bmi_cat"] = pd.Series(np.select([df.bmi < 25, df.bmi < 30], ["normal", "over"], default="obese"), index=df.index).where(df.bmi.notna())
    # smoking as coded in Table 1: >= 100 cigarettes lifetime yes / no
    df["smoke100"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    df["alc_days"] = _alcohol_days(lib, df)
    # hypertension: told, mean SBP > 130 or DBP > 80 ("exceeding"), or taking prescription
    pos = (df.BPQ020 == 1) | (df.BPQ_med_now == 1) | (df.sbp > 130) | (df.dbp > 80)
    known = df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna())
    df["htn"] = pd.Series(np.where(pos, 1.0, np.where(known, 0.0, np.nan)), index=df.index)
    # diabetes: self-report, insulin or other diabetes treatment, HbA1c > 6.5, fasting glucose >= 7.0 mmol/L (126 mg/dL)
    pos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (df.DIQ070 == 1) | (df.hba1c > 6.5) | (df.glucose >= 126)
    known = df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna()
    df["dm"] = pd.Series(np.where(pos, 1.0, np.where(known, 0.0, np.nan)), index=df.index)
    terms = ["x", "female", "age"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "bmi_cat", "normal")
    terms += ["waist"]
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["smoke100", "alc_days", "sbp", "dbp", "ast", "alt", "uric_acid", "tc", "ldl", "hdl", "calcium", "htn", "dm"]
    df["w"] = lib.pooled_weight(df, "fasting")
    cancer = lib.raw(df, lib.stems("MCQ"), "MCQ220") == 1
    dom = (df.age >= 20) & ~cancer & (df.pregnant != 1) & df.uacr.notna() & e.notna() & df.x.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "cycles (spec ambiguity)": "reading_1: 1999-2000 to 2015-2016 plus 2017-March 2020 pre-pandemic file (P_)",
                "eGFR equation (spec ambiguity)": "reading_1: CKD-EPI 2009 with race term",
                "inflection point (spec ambiguity)": "reading_1, not used by M1",
                "weight": "unnamed survey weight; fasting subsample weight (TyG needs fasting glucose and TG); WTSAF4YR for 1999-2002, WTSAFPRP for 2017-2020",
                "age": "'< 20 years of age' excluded read as age >= 20 kept",
                "cancer": "MCQ220 = 1 excluded; missing MCQ220 kept",
                "missing covariates": "the paper used mode/median imputation; prereg section 3 fixes complete-case, so no imputation",
                "age coding": "continuous (text reports mean age)",
                "BMI coding": "Table 1 categories normal < 25, overweight 25 to < 30, obese >= 30 (the only coding the paper reports)",
                "WC": "continuous cm",
                "race": "RIDRETH1 5 groups as in Table 1",
                "education": "Table 1: less than high school / high school or GED / above high school (educ3); 'Others' (7/9) missing",
                "smoking": "Table 1 coding: >= 100 cigarettes lifetime (SMQ020) yes/no",
                "alcohol": "days of alcohol use in past year (ALQ120Q x unit); lifetime non-drinkers 0; 2017-2020 ALQ121 categories mapped to days/year (every day 365, nearly every day 312, 3-4/wk 182, 2/wk 104, 1/wk 52, 2-3/mo 30, 1/mo 12, 7-11/yr 9, 3-6/yr 4.5, 1-2/yr 1.5, none 0)",
                "hypertension": "BPQ020 = 1 or BPQ050A = 1 (taking prescription) or mean SBP > 130 or DBP > 80 (literal 'exceeding')",
                "diabetes": "DIQ010 = 1 or DIQ050 = 1 or DIQ070 = 1 or HbA1c > 6.5 (literal '>') or fasting glucose >= 126 mg/dL",
                "lab covariates": "AST, ALT, uric acid, TC, LDL-C (fasting file), HDL-C, total calcium (BIOPRO) continuous in reported units",
            },
            "notes": "written before first execution; no post-run changes"}

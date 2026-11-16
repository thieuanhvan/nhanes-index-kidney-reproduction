"""M1 builder, 69700c2b09 (Lipids Health Dis 2024, doi [DOI withheld]). NHHR continuous, DKD among
adults with diabetes, single adjusted model (Table 2, OR 1.07).
Spec: NHANES 2005-2016, age >= 18 with diabetes, not pregnant, DKD = UACR >= 30 mg/g or eGFR < 60, weight WTDR2D / 6.
Sources (fulltext_txt/69700c2b09.txt): exclusions line 13, NHHR line 16 (formula image M1 in the XML: non-HDL-C / HDL-C),
DKD line 19, covariates line 22, weights line 25, Table 1 coding lines 120-172, Table 2 footnote line 186.
"""
import numpy as np, pandas as pd

REC_ID = "69700c2b09"
CYCLES = ["2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016"]


def _yn(s):
    return (s == 1).astype(float).where(s.isin([1, 2]))


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["NHHR"][0](df)
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())
    # weight: dietary two-day sample weight / 6
    d1 = {c: "DR1TOT" + s for c, s in lib.SUF.items()}
    d2 = {c: "DR2TOT" + s for c, s in lib.SUF.items()}
    df["w"] = lib.raw(df, d1, "WTDR2D").fillna(0) / 6
    # marital status (Table 1): married or living with partner vs widowed/divorced/separated/never married
    mar = lib.raw(df, lib.stems("DEMO"), "DMDMARTL")
    df["married"] = mar.isin([1, 6]).astype(float).where(mar.isin([1, 2, 3, 4, 5, 6]))
    # family PIR (Table 1 categories): <= 1.0, 1.0 to 3.0, > 3.0
    df["pir_cat"] = pd.Series(np.select([df.pir <= 1, df.pir <= 3], ["le1", "1to3"], default="gt3"), index=df.index).where(df.pir.notna())
    # BMI (Table 1 categories): < 25, 25 to 30, >= 30
    df["bmi_cat"] = pd.Series(np.select([df.bmi < 25, df.bmi < 30], ["lt25", "25to30"], default="ge30"), index=df.index).where(df.bmi.notna())
    # smoking: > 100 cigarettes in life (SMQ020); drinking: >= 12 drinks in any one year (ALQ101)
    df["smoke100"] = _yn(df.SMQ020)
    df["drinker"] = _yn(lib.raw(df, lib.stems("ALQ"), "ALQ101"))
    # history of hypertension: self-report or BP exceeding 140/90
    pos = (df.BPQ020 == 1) | (df.sbp > 140) | (df.dbp > 90)
    known = df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna())
    df["htn"] = pd.Series(np.where(pos, 1.0, np.where(known, 0.0, np.nan)), index=df.index)
    # hyperuricemia: SUA > 420 umol/L men, > 360 umol/L women (mg/dL x 59.48)
    ua = df.uric_acid * 59.48
    df["hyperuric"] = pd.Series(np.where(df.female == 1, ua > 360, ua > 420), index=df.index).astype(float).where(ua.notna() & df.female.notna())
    # physical activity: vigorous leisure-time activity >= 10 min (PAD200 in 2005-2006, PAQ650 from 2007)
    pa = lib.raw(df, lib.stems("PAQ"), ["PAQ650", "PAD200"])
    df["active"] = (pa == 1).astype(float).where(pa.isin([1, 2, 3]))
    # dietary cholesterol: mean of the two recalls (mg)
    c1 = lib.raw(df, d1, "DR1TCHOL")
    c2 = lib.raw(df, d2, "DR2TCHOL")
    df["diet_chol"] = (c1 + c2) / 2
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["married"]
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "pir_cat", "le1")
    terms += lib.dummies(df, "bmi_cat", "lt25")
    terms += ["smoke100", "drinker", "htn", "hyperuric", "active", "diet_chol", "ast", "alt"]
    dom = ((df.age >= 18) & (df.diabetes_any == 1) & (df.pregnant != 1) & df.x.notna() & e.notna() & df.uacr.notna()
           & (df.w > 0))
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "weight (spec ambiguity)": "reading_1: WTDR2D / 6 as stated",
                "NHHR": "(TC - HDL-C) / HDL-C (formula image in the XML; standard definition)",
                "exposure scale": "continuous NHHR per 1 unit",
                "age": "'under the age of 18' excluded read as age >= 18 kept",
                "diabetes (population)": "DIQ010 = 1 or DIQ050/DIQ070 = 1 or fasting glucose >= 126 or HbA1c >= 6.5 (diabetes_any; borderline = no)",
                "missing weighted data": "WTDR2D missing or 0 excluded (domain requires w > 0)",
                "eGFR equation": "CKD-EPI 2009 with race (not stated)",
                "age coding": "continuous (Table 1 median)",
                "race": "RIDRETH1 5 groups as in Table 1",
                "marital status": "Table 1: married or living with partner (DMDMARTL 1, 6) vs widowed/divorced/separated/never married; refused missing",
                "education": "Table 1: less than high school / high school / college or above (educ3, DMDEDUC2, asked at 20+; 18-19 missing)",
                "family PIR": "Table 1 categories <= 1.0, 1.0 to 3.0, > 3.0",
                "BMI": "Table 1 categories < 25, 25 to 30, >= 30",
                "smoking": "SMQ020 >= 100 cigarettes yes/no",
                "alcohol": "ALQ101 >= 12 drinks in any one year yes/no ('at least 12 ... in any given year')",
                "hypertension": "BPQ020 = 1 or mean SBP > 140 or DBP > 90 (literal 'exceeding'); medication not part of the stated definition",
                "hyperuricemia": "SUA > 420 umol/L men, > 360 umol/L women ('mmol/L' in text read as umol/L)",
                "physical activity": "vigorous leisure activity: PAD200 = 1 (2005-2006), PAQ650 = 1 (2007-2016); else inactive; 7/9 missing",
                "dietary cholesterol": "mean of DR1TCHOL and DR2TCHOL, mg; both days required (WTDR2D > 0 implies two valid recalls)",
                "AST, ALT": "continuous U/L",
            },
            "notes": "written before first execution; no post-run changes"}

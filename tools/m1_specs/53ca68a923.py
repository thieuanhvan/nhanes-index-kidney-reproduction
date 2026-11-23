"""M1 builder, 53ca68a923 (CALLY index and CKD; source withheld).
Spec: CALLY Q4 vs Q1, CKD = eGFR < 60 and/or UACR >= 30, Model 3, NHANES 1999-2000 to 2009-2010, survey weighted.
Sources in fulltext_txt/53ca68a923.txt: population line 12, index and CKD line 15, covariates line 21, models and
weighting line 24, Table I lines 68-100 (Smokers yes/no, Married, Above high school), Table III.
"""
REC_ID = "53ca68a923"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010"]


def build(lib):
    df = lib.table(CYCLES)
    df["cally"] = lib.INDEX["CALLY_gL"][0](df)
    df["y"] = ((df.egfr_ckdepi2009 < 60) | (df.uacr >= 30)).astype(float).where(df.uacr.notna() & df.egfr_ckdepi2009.notna())
    s = lib.stems
    mar = lib.raw(df, s("DEMO"), "DMDMARTL")
    df["married"] = (mar == 1).astype(float).where(mar.isin([1, 2, 3, 4, 5, 6]))
    df["educ_above_hs"] = df.educ_raw.isin([4, 5]).astype(float).where(df.educ_raw.notna())
    df["ever_smoker"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    pills = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    dpos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (pills == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5)
    dobs = df.DIQ010.isin([1, 2, 3]) | df.glucose.notna() | df.hba1c.notna()
    df["diab"] = dpos.astype(float).where(dpos | dobs)
    dom = (df.age >= 20) & (df.pregnant != 1) & df.cally.notna() & (df.cally > 0) & df.y.notna()
    qn, cuts = lib.quantile_dummies(df, "cally", dom, k=4, prefix="cally_q")
    terms = qn + ["age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["married", "pir", "educ_above_hs", "ever_smoker", "hypertension_any_140", "diab", "cvd_history",
              "bmi", "tg", "hdl", "ldl", "uric_acid"]
    df["w"] = lib.pooled_weight(df, "fasting")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "cally_q4", "weight": "w",
            "readings": {
                "estimate": "Q4 vs Q1 Model 3 (reading_1); quartiles of CALLY cut unweighted on the domain (ln does not change cuts)",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated, race mentioned)",
                "UACR threshold": ">= 30 mg/g as written",
                "race": "RIDRETH1 5 groups",
                "marital": "married (DMDMARTL 1) vs unmarried (2-6); 77/99 missing",
                "education": "above high school (DMDEDUC2 4-5) vs not",
                "smoking history": "binary ever smoker (SMQ020, 100 cigarettes): Table I reports Smokers yes/no and Methods says former/current users",
                "hypertension": "hypertension_any_140: self-report, mean BP >= 140/90, BPQ050A medication",
                "diabetes": "DIQ010 == 1, FPG >= 126, HbA1c >= 6.5, insulin DIQ050, pills DIQ070/DID070",
                "CVDs": "cvd_history MCQ160B-F (CHF, CHD, angina, MI, stroke)",
                "continuous covariates": "PIR, BMI, TG, HDL, LDL, uric acid continuous (TC listed in Methods but not in the Model 3 list, not used)",
                "weight": "weighting variable not named; pooled fasting-subsample weight because TG and LDL come from the fasting subsample (brief rule 8)"},
            "notes": f"builder by agent 15/11; quartile cuts on domain {[int(c * 10000) / 10000 for c in cuts]}"}

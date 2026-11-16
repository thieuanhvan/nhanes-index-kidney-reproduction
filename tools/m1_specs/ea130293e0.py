"""M1 builder, ea130293e0 (RA = RDW/albumin and DKD, Sci Rep 2026, doi [DOI withheld]).
Spec: RA continuous per 1 dL/g, DKD = eGFR < 60 or UACR >= 30 in diabetes, NHANES 2011-2018, Model 3 (reading_1),
survey weighted. Sources in fulltext_txt/ea130293e0.txt: population (Study design, line 13), diabetes and DKD
(Endpoint definition, line 21), covariates (line 25), Model 3 (Statistical analysis, line 28), Scr and eGFR in the
fully adjusted model (multicollinearity paragraph).
"""
import numpy as np

REC_ID = "ea130293e0"
CYCLES = ["2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["RAR"][0](df)
    # diabetes: HbA1c >= 6.5, FPG >= 7.0 mmol/L (126 mg/dL), 2-h OGTT >= 11.1 mmol/L (200 mg/dL), self-reported
    # diagnosis, insulin or oral hypoglycaemic agents
    ogtt = lib.raw(df, s("OGTT"), "LBXGLT")
    diq070 = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    pos = ((df.DIQ010 == 1) | (df.DIQ050 == 1) | (diq070 == 1) | (df.hba1c >= 6.5) | (df.glucose >= 126) | (ogtt >= 200))
    df["diab"] = pos.astype(float)
    # outcome: eGFR < 60 (CKD-EPI, version not stated -> 2009 with race) or UACR >= 30
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())
    # hypertension: physician diagnosis or measured SBP >= 140 and/or DBP >= 90
    known = df.dx_hypertension.notna() | (df.sbp.notna() & df.dbp.notna())
    df["htn"] = ((df.dx_hypertension == 1) | (df.sbp >= 140) | (df.dbp >= 90)).astype(float).where(known)
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "smoking", "never")
    terms += ["bmi", "htn", "hdl", "ldl", "scr", "bun", "uric_acid", "egfr_ckdepi2009"]
    df["w"] = lib.pooled_weight(df, "fasting")
    dom = ((df.age >= 18) & (df.diab == 1) & df.rdw.notna() & df.albumin_serum.notna()
           & df.uacr.notna() & df.egfr_ckdepi2009.notna())
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "Model 3 covariate list": "reading_1: Model 2 (age, sex, race, education, smoking) plus BMI, blood pressure, lipid profile, renal function parameters",
                "main estimate cohort": "reading_1: NHANES",
                "inflection point": "reading_1: none",
                "age": "excluding those under 18 read as age >= 18",
                "diabetes": "paper composite: DIQ010 = 1, insulin DIQ050, pills DIQ070, HbA1c >= 6.5, FPG >= 126 mg/dL, OGTT >= 200 mg/dL (OGTT 2011-2016 only, not run in 2017-2018)",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated)",
                "UACR threshold": ">= 30 mg/g as written",
                "race": "RIDRETH1 5 groups (coding not stated)",
                "education": "educ3 (less than high school, high school, more than high school; coding not stated)",
                "smoking": "never, former, current as defined in Covariates",
                "BMI": "continuous (Methods give no categories; Table 1 categories are descriptive)",
                "blood pressure": "hypertension as defined in Covariates: diagnosis or SBP >= 140 or DBP >= 90",
                "lipid profile": "HDL-C and LDL-C (the lipids named in Covariates)",
                "renal function parameters": "Scr, BUN, serum uric acid (named in Covariates) and eGFR (named as in the fully adjusted model in the multicollinearity paragraph)",
                "etc. in Model 3": "not extended beyond the named clinical groups; PIR not added (socioeconomic, not a clinical indicator, absent from Model 2)",
                "weight": "weighted; fasting subsample weight because LDL-C (fasting subsample) is a model term",
                "missing data": "paper used multiple imputation; M1 is complete case by design",
            },
            "notes": "LDL-C restricts the complete cases to the fasting subsample, hence fasting weights (brief rule 8)"}

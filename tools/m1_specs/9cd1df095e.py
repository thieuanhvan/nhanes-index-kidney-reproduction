"""M1 builder, 9cd1df095e (source withheld). NHR per 1 unit, albuminuria, Model 3.
Spec: NHANES 2017-March 2020 pre-pandemic file, adults >= 20, albuminuria ACR > 30 mg/g, survey weights (not named),
Table 2 Model 3 OR 1.13.
Sources (fulltext_txt/9cd1df095e.txt): population and exclusions (Methods, Study population), exposure and outcome
(Exposure variable and outcome definition), covariates paragraph, Table 1 coding, Table 2 footnote (Model 3 list).
"""
import numpy as np, pandas as pd

REC_ID = "9cd1df095e"
CYCLES = ["2017-2020"]


def _yn(s):
    return (s == 1).astype(float).where(s.isin([1, 2]))


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["NHR"][0](df)
    df["y"] = (df.uacr > 30).astype(float).where(df.uacr.notna())
    st = lib.stems
    # marital status as in Table 1 ("Marriage yes/no"): DMDMARTZ 1 married/living with partner vs 2, 3
    mar = lib.raw(df, {"2017-2020": "P_DEMO"}, "DMDMARTZ")
    df["married"] = (mar == 1).astype(float).where(mar.isin([1, 2, 3]))
    # physical activity status as in Table 1 ("Moderate work activity yes/no"): PAQ620
    df["pa_modwork"] = _yn(lib.raw(df, st("PAQ"), "PAQ620"))
    # alcohol consumption yes/no: ALQ111 ever had a drink of any kind of alcohol
    df["alcohol"] = _yn(lib.raw(df, st("ALQ"), "ALQ111"))
    # smoking status yes/no: SMQ020 >= 100 cigarettes in life
    df["smoke100"] = _yn(df.SMQ020)
    # history of hyperlipidemia: BPQ080 told high cholesterol; kidney stones: KIQ026
    df["hyperlip"] = _yn(lib.raw(df, st("BPQ"), "BPQ080"))
    df["kstone"] = _yn(lib.raw(df, st("KIQ_U"), "KIQ026"))
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["pir"]
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["pa_modwork", "married", "alcohol", "smoke100", "sbp", "dbp", "scr", "bun", "uric_acid", "glucose",
              "hba1c", "tg", "ldl", "egfr_ckdepi2009", "bmi", "dx_diabetes", "dx_hypertension", "hyperlip", "kstone"]
    df["w"] = lib.pooled_weight(df, "fasting")
    dom = (df.age >= 20) & df.x.notna() & df.uacr.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "Model 3 covariate list (spec ambiguity)": "reading_1: Table 2 footnote list, includes serum creatinine",
                "weight and pooling (spec ambiguity)": "reading_1: unnamed survey weights, no pooling rule; the smallest subsample used (fasting glucose, TG, LDL-C are covariates) gives the fasting weight WTSAFPRP",
                "age": "'aged < 20 years' excluded read as age >= 20 kept",
                "missing covariates": "the paper used median/mode imputation; prereg fixes complete-case, so no imputation",
                "NHR": "neutrophil count (10^3/uL) / HDL-C (mmol/L, mg/dL / 38.67)",
                "race": "RIDRETH1 5 groups as in Table 1 (Other race includes Asian)",
                "education": "Table 1: less than high school / high school / more than high school (educ3)",
                "marital status": "Table 1 'Marriage yes/no': DMDMARTZ = 1 (married or living with partner) yes, 2/3 no, 77/99 missing",
                "physical activity status": "Table 1 row 'Moderate work activity yes/no': PAQ620",
                "alcohol consumption": "yes/no, not defined in the text; ALQ111 ever had a drink of alcohol",
                "smoking status": "yes/no; SMQ020 >= 100 cigarettes lifetime",
                "diabetes, hypertension": "'history of': self-reported diagnosis DIQ010 = 1 (borderline = no), BPQ020 = 1",
                "hyperlipidemia": "'history of': BPQ080 told high blood cholesterol",
                "kidney stones": "KIQ026 ever had kidney stones",
                "glomerular filtration rate": "CKD-EPI 2009 with race (equation not stated)",
                "glucose, TG, LDL-C": "fasting plasma glucose LBXGLU, fasting TG LBXTR and LDL-C LBDLDL (P_GLU, P_TRIGLY)",
                "continuous covariates": "age, PIR, BMI, SBP, DBP, Scr, BUN, uric acid, glucose, HbA1c, TG, LDL-C, eGFR as continuous (Table 1 reports them as means)",
            },
            "notes": "written before first execution; no post-run changes"}

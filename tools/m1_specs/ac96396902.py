"""M1 builder, ac96396902 (ALBI and CVD-CKD comorbidity; source withheld).
Spec (reading_1): CKD alone (UACR >= 30 or eGFR < 60), ALBI continuous per unit, adjusted sensitivity model of the
Results text (Supplementary Table 2), NHANES 2011-2018, sample weights.
Sources in fulltext_txt/ac96396902.txt: participants (line 16), ALBI (line 19), CVD and CKD (line 22),
covariates (line 25), weights and models (line 28), CKD-alone adjustment list (Sensitivity analysis, line 45),
Table 1 (education levels, smoking, drinking, hypertension coding).
"""
import numpy as np
import pandas as pd

REC_ID = "ac96396902"
CYCLES = ["2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["ALBI"][0](df)
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())
    # nondrinker: answered no to at least 12 drinks in life or in any one year (2011-2016 ALQ101 and ALQ110 both no);
    # 2017-2018 has neither item: ALQ111 (ever had a drink) used as the lifetime item
    a = s("ALQ")
    a101, a110, a111 = lib.raw(df, a, "ALQ101"), lib.raw(df, a, "ALQ110"), lib.raw(df, a, "ALQ111")
    yes = (a101 == 1) | (a110 == 1) | (a111 == 1)
    no = ((a101 == 2) & (a110 == 2)) | (a111 == 2)
    df["drinker"] = pd.Series(np.where(yes, 1.0, np.where(no, 0.0, np.nan)), index=df.index)
    # smoker: at least 100 cigarettes in lifetime
    df["smoker"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    # education: Table 1 five levels (DMDEDUC2 1-5)
    df["educ5"] = df.educ_raw.where(df.educ_raw.isin([1, 2, 3, 4, 5]))
    terms = ["x", "female", "age"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ5", 1.0)
    terms += ["pir", "bmi", "drinker", "smoker", "dx_hypertension"]
    df["w"] = lib.pooled_weight(df, "mec")
    dom = (df.age >= 20) & df.cvd_history.notna() & df.uacr.notna() & e.notna() & df.x.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "outcome of the main estimate": "reading_1: CKD alone, OR 2.10 per unit (Results text, Supplementary Table 2)",
                "diabetes in the adjustment set": "reading_1: not included",
                "inflection point": "reading_1: none",
                "age": "adults read as age >= 20 (the 5 education levels are DMDEDUC2, asked of age 20 and over)",
                "CVD information complete": "CHD, CHF, MI, angina, stroke items all answered (cvd_history not missing)",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated)",
                "UACR threshold": ">= 30 mg/g as written",
                "race": "5 groups of Table 1 (RIDRETH1)",
                "education": "5 levels of Table 1 (DMDEDUC2), reference less than 9th grade",
                "PIR, BMI, age": "continuous (Table 1 means)",
                "drinking": "nondrinker = no to 12 drinks in any one year and no to 12 drinks in lifetime (2011-2016); 2017-2018 ALQ111 ever drank",
                "smoking": "smoker = at least 100 cigarettes in lifetime (SMQ020 = 1), as written; Table 1 shows 19 % yes, which fits current smoking, but the text is followed",
                "hypertension": "self-reported diagnosis (BPQ020), as written",
                "weight": "MEC weight pooled over 4 cycles (variable not named)",
            },
            "notes": "ALBI per INDEX['ALBI']: log10 bilirubin umol/L * 0.66 - 0.085 * albumin g/L"}

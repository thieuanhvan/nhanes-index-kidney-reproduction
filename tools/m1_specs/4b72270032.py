"""M1 builder, 4b72270032 (RAR and CKD; source withheld).
Spec: RAR continuous per unit, CKD = eGFR < 60 or UACR >= 30, NHANES 2005-2018, Table 2 Model 3 footnote, weighted.
Sources in fulltext_txt/4b72270032.txt: population (2.1, line 14), CKD (2.2, line 17), covariate definitions (2.4,
lines 23-25), weights (2.5, line 28), Model 3 list (Table 2 footnote, line 153).
"""
import numpy as np
import pandas as pd

REC_ID = "4b72270032"
CYCLES = ["2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["RAR"][0](df)
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())
    # diabetes: physician diagnosis, HbA1c > 6.5 %, random blood glucose >= 11.1 mmol/L (200 mg/dL, serum glucose of
    # the biochemistry profile, not fasting), diabetes medication or insulin
    diq070 = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    pos = ((df.DIQ010 == 1) | (df.DIQ050 == 1) | (diq070 == 1) | (df.hba1c > 6.5) | (df.glucose_serum >= 200))
    known = df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose_serum.notna()
    df["diab"] = pos.astype(float).where(known)
    # PIR in 4 groups: <= 1.30, 1.31-1.85, 1.86-3.50, > 3.50
    p = df.pir
    df["pir4"] = pd.Series(np.select([p <= 1.30, p <= 1.85, p <= 3.50, p > 3.50], ["a_le130", "b_131_185", "c_186_350", "d_gt350"], default=""),
                           index=df.index).where(p.notna())
    terms = ["x", "age", "female", "bmi"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "pir4", "a_le130")
    terms += lib.dummies(df, "smoking", "never")
    terms += ["hypertension_any_130", "cvd_history", "diab", "bun", "alt", "ast", "hba1c", "glucose_serum"]
    df["w"] = lib.pooled_weight(df, "mec")
    dom = (df.age >= 20) & df.scr.notna() & df.uacr.notna() & df.x.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "abstract has no numeric estimate": "reading_1: continuous RAR, Table 2 Model 3, OR 1.58",
                "Model 3 covariate set": "reading_1: Table 2 footnote list (with HbA1c and glucose; no alcohol, no lipids)",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated)",
                "UACR threshold": ">= 30 mg/g as written",
                "race": "5 groups as written (RIDRETH1)",
                "education": "3 groups as written (educ3); 'no record' rows treated as missing",
                "PIR": "4 groups as written in 2.4 (<= 1.30, 1.31-1.85, 1.86-3.50, > 3.50)",
                "age, BMI": "continuous (Methods give no categories; Table 1 groups are descriptive)",
                "smoking": "never, former, current as written",
                "hypertension": "diagnosis, mean BP >= 130/80 or current antihypertensive medication (hypertension_any_130)",
                "diabetes": "diagnosis, HbA1c > 6.5 (strict as written), random serum glucose >= 200 mg/dL, DIQ070/DID070 pills or DIQ050 insulin",
                "CVD": "CHF, CHD, angina, MI or stroke (MCQ160B-F)",
                "glucose covariate": "serum glucose of the biochemistry profile (Table 1 glucose is reported for the whole sample; fasting glucose would cover only the fasting subsample)",
                "weight": "MEC weight pooled over 7 cycles (no fasting-only term)",
            },
            "notes": "units of BUN, glucose kept in NHANES conventional units (scale does not change the RAR OR)"}

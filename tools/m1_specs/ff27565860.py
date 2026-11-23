"""M1 builder, ff27565860 (NHHR and DKD in diabetes; source withheld).
Spec: NHHR continuous per 1 unit, Table 2 Model 3, DKD = diabetes plus UACR >= 30 and/or eGFR < 60, NHANES 1999-2020
(2017-March 2020 pre-pandemic release used for the last period).
Sources (fulltext_txt/ff27565860.txt): population (criteria only in Fig 1), exposure, outcome and diabetes definition,
covariates paragraph (smoking, BMI and PIR categories), statistics paragraph, Table 1 coding, Table 2.
"""
import numpy as np
import pandas as pd

REC_ID = "ff27565860"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010",
          "2011-2012", "2013-2014", "2015-2016", "2017-2020"]


def flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def labels(index, pairs):
    out = pd.Series(np.nan, index=index, dtype=object)
    for m, lab in reversed(pairs):
        out[m.fillna(False).to_numpy(bool)] = lab
    return out


def yn(lib, df, base, col):
    v = lib.raw(df, lib.stems(base), col)
    return flag(v == 1, v.isin([1, 2]))


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["NHHR"][0](df)
    pills = lib.raw(df, lib.stems("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, lib.stems("DIQ"), "DIQ050")
    diab = (df.dx_diabetes == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5) | (ins == 1) | (pills == 1)
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())
    dom = (df.age >= 20) & diab & df.x.notna()
    # marital status, 6 categories (Table 1); P_DEMO has only the 3-group DMDMARTZ, so the 2017-2018 part of the
    # 2017-2020 release takes DMDMARTL from DEMO_J (same SEQN); 2019-March 2020 participants stay missing
    mar = lib.raw(df, lib.stems("DEMO"), "DMDMARTL")
    pj = (df.cycle == "2017-2020").to_numpy()
    dj = lib._read("2017-2018", "DEMO_J")
    if dj is not None and pj.any():
        mar[pj] = dj.drop_duplicates("SEQN").set_index("SEQN")["DMDMARTL"].reindex(df.loc[pj, "SEQN"]).to_numpy(float)
    df["marital6"] = mar.map({1: "married", 2: "widowed", 3: "divorced", 4: "separated", 5: "never", 6: "partner"})
    df["educ5"] = df.educ_raw.map({1: "lt9", 2: "9to11", 3: "hs", 4: "somecoll", 5: "coll"})
    df["bmi3"] = labels(df.index, [(df.bmi <= 25, "le25"), (df.bmi <= 30, "25to30"), (df.bmi > 30, "gt30")])
    df["pir3"] = labels(df.index, [(df.pir <= 1, "le1"), (df.pir <= 3, "1to3"), (df.pir > 3, "gt3")])
    df["htn_dx"] = df["dx_hypertension"]
    df["stroke"] = yn(lib, df, "MCQ", "MCQ160F")
    df["cancer"] = yn(lib, df, "MCQ", "MCQ220")
    df["chd"] = yn(lib, df, "MCQ", "MCQ160C")
    df["hichol"] = yn(lib, df, "BPQ", "BPQ080")
    df["liver"] = yn(lib, df, "MCQ", "MCQ160L")
    terms = ["x"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["female", "age"]
    terms += lib.dummies(df, "marital6", "married")
    terms += lib.dummies(df, "educ5", "lt9")
    terms += lib.dummies(df, "bmi3", "le25")
    terms += lib.dummies(df, "pir3", "le1")
    terms += ["htn_dx", "stroke", "cancer", "chd", "hichol", "liver"]
    terms += lib.dummies(df, "smoking", "never")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": None,
            "readings": {
                "estimate (spec ambiguity 1)": "reading_1: Table 2 Model 3 continuous NHHR OR 1.06",
                "Model 3 covariates (spec ambiguity 2)": "reading_1: covariates paragraph plus smoking; HDL-C and TC not added",
                "survey weights (spec ambiguity 3)": "reading_1: unweighted",
                "cycles": "1999-2000 to 2015-2016 plus the 2017-March 2020 pre-pandemic release (P_ files)",
                "population": "Fig 1 criteria not in the text: adults age >= 20 with diabetes as defined in the paper; pregnancy not excluded (not stated); missing NHHR, UACR, eGFR left to complete case",
                "diabetes": "DIQ010 = 1, fasting glucose >= 126 mg/dL (7.0 mmol/L), HbA1c >= 6.5, insulin (DIQ050) or pills (DIQ070; DID070 in 2005-2008)",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated in the spec; the printed formula is the 2009 equation without its race rows)",
                "UACR threshold": ">= 30 mg/g",
                "race": "RIDRETH1 5 groups (Table 1), reference NH white",
                "age": "continuous (Table 1 mean +- SD)",
                "marital status": "DMDMARTL 6 categories (Table 1), reference married; 2017-2020 from DEMO_J where the SEQN is a 2017-2018 participant, otherwise missing (P_DEMO has only 3 groups)",
                "education": "DMDEDUC2 5 categories (Table 1), reference less than 9th grade",
                "BMI": "3 categories as written: <= 25, 25 to 30, > 30 (methods state the categories; Table 1 shows the mean)",
                "PIR": "3 categories as written: <= 1, 1 to 3, > 3",
                "hypertension": "questionnaire diagnosis BPQ020 = 1 ('medical history ... from NHANES questionnaires')",
                "stroke, cancer, coronary heart disease, liver disease": "MCQ160F, MCQ220, MCQ160C, MCQ160L = 1",
                "hypercholesterolemia": "told high cholesterol BPQ080 = 1",
                "smoking": "never / former / current (SMQ020, SMQ040), reference never"},
            "notes": "2017-2020 release used instead of 2017-2018 (spec cycles say 2017-2020 pre-pandemic)"}

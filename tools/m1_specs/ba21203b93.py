"""M1 builder, ba21203b93 (PSA indicators and DKD in men with diabetes, Eur J Med Res 2025, doi [DOI withheld]).
Spec: FPSA/TPSA ratio (x 10^2) tertile T3 vs T2 (reference T2), CKD in men with diabetes, Model IV (15 covariates),
NHANES 2001-2010, no weighting.
Sources in the paper (fulltext_txt/ba21203b93.txt): population line 17, DM and CKD line 21, covariate coding line 24,
models line 28, Table 1 categories lines 96-136, Table 2 reference levels lines 140-182, Table 3 footnote line 205.
"""
import numpy as np
import pandas as pd

REC_ID = "ba21203b93"
CYCLES = ["2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010"]
OGTT = {"2005-2006": "OGTT_D", "2007-2008": "OGTT_E", "2009-2010": "OGTT_F"}


def _miss(v, codes=(7, 9, 77, 99)):
    return v.where(~v.isin(codes))


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["FPSA_TPSA"][0](df).replace([np.inf, -np.inf], np.nan)
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())

    # DM: HbA1c >= 6.5; OGTT >= 200, FPG >= 126 or random glucose >= 200; diagnosis; insulin or other antidiabetic medication
    ogtt = lib.raw(df, OGTT, "LBXGLT")
    df["dm"] = ((df.hba1c >= 6.5) | (ogtt >= 200) | (df.glucose >= 126) | (df.glucose_serum >= 200) | (df.DIQ010 == 1)
                | (df.DIQ050 == 1) | (df.DIQ070 == 1)).astype(float)
    # prostate cancer: MCQ230A-D cancer type code 30 (prostate)
    pc = pd.Series(False, index=df.index)
    for c in ["MCQ230A", "MCQ230B", "MCQ230C", "MCQ230D"]:
        pc |= lib.raw(df, s("MCQ"), c) == 30

    df["age_gt60"] = (df.age > 60).astype(float).where(df.age.notna())
    df["educ_gt_hs"] = (df.educ3 == "gt_hs").astype(float).where(df.educ3.notna())
    mar = _miss(lib.raw(df, s("DEMO"), "DMDMARTL"))
    df["married"] = mar.isin([1, 6]).astype(float).where(mar.isin([1, 2, 3, 4, 5, 6]))
    df["bmi_cat"] = pd.Series(np.select([df.bmi < 25, df.bmi <= 30, df.bmi > 30], ["lt25", "25-30", "gt30"], default=""),
                              index=df.index).astype(object).where(df.bmi.notna(), np.nan)
    # smokers: >= 100 cigarettes in lifetime
    df["smoke_yes"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    # drinkers: 12 or more drinks per year (ALD100 2001-2002, ALQ101 2003-2010)
    a = _miss(lib.raw(df, s("ALQ"), ["ALQ101", "ALD100"]))
    df["drink_yes"] = (a == 1).astype(float).where(a.isin([1, 2]))
    df["stroke"] = (df.MCQ160F == 1).astype(float).where(df.MCQ160F.isin([1, 2]))
    df["chd"] = (df.MCQ160C == 1).astype(float).where(df.MCQ160C.isin([1, 2]))
    pos = (df.sbp >= 140) | (df.dbp >= 90) | (df.BPQ020 == 1)
    df["htn"] = pd.Series(np.where(pos, 1.0, np.where(df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna()), 0.0, np.nan)), index=df.index)

    dom = ((df.female == 0) & (df.dm == 1) & (df.age >= 20) & ~pc & df.x.notna() & df.uacr.notna() & e.notna())
    _, cuts = lib.quantile_dummies(df, "x", dom, k=3)
    g = df["x_qgroup"]
    df["x_t1"] = (g == 1).astype(float).where(g.notna())
    df["x_t3"] = (g == 3).astype(float).where(g.notna())

    terms = ["x_t1", "x_t3", "age_gt60"]
    terms += lib.dummies(df, "race1", "mexican_american")
    terms += ["educ_gt_hs", "married"]
    terms += lib.dummies(df, "bmi_cat", "lt25")
    terms += ["waist", "smoke_yes", "drink_yes", "stroke", "chd", "htn", "tc", "tg_serum", "hdl", "hba1c"]

    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x_t3", "weight": None,
            "readings": {
                "main estimate": "reading_1 FPSA/TPSA ratio, Model IV, T3 vs T2 OR 2.16",
                "covariate list": "reading_1 Table 3 footnote Model IV (15 covariates)",
                "continuous vs category": "reading_1 tertile contrast T3 vs T2",
                "tertiles": "unweighted tertiles of FPSA/TPSA x 100 cut on the domain (printed ranges T1 4-21, T2 24-32, T3 36-51 not used as fixed cuts)",
                "eGFR equation": "not stated: CKD-EPI 2009 with race",
                "CKD": "UACR >= 30 and/or eGFR < 60 in men with diabetes",
                "survey weighting": "none stated: unweighted",
                "age limit": "none stated (PSA measured in men >= 40); age >= 20 kept from recon CFG, no effect",
                "prostate cancer": "MCQ230A-D = 30",
                "diabetes": "as written; random glucose = serum LBXSGL >= 200; OGTT 2005-2010 only; pills DIQ070 absent 2005-2008",
                "age": "binary <= 60 / > 60 (Tables 1-2)",
                "race": "RIDRETH1 5 groups, reference Mexican American (Table 2)",
                "education": "> high school (DMDEDUC2 4-5) vs <= high school (Methods)",
                "marital status": "binary married or living with partner vs never married or others (Tables 1-2)",
                "BMI": "3 groups < 25, 25-30, > 30",
                "smoking": "smoker = >= 100 cigarettes in lifetime (SMQ020 = 1), as defined in Methods",
                "drinking": "drinker = 12 or more drinks in a year (ALD100 / ALQ101 = 1)",
                "heart disease": "coronary heart disease MCQ160C (Methods 'history of coronary heart disease')",
                "hypertension": "SBP >= 140 and/or DBP >= 90, or physician diagnosis (BPQ020)",
                "TG": "serum biochemistry triglyceride LBXSTR: only 4 participants were excluded for missing TG, which rules out the fasting-subsample LBXTR",
                "TC, HDL, HbA1c, waist": "continuous",
                "missing data": "paper used mean/median imputation; M1 fit is complete case (automatic)",
            },
            "notes": "builder written before first execution"}

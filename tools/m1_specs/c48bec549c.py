"""M1 builder, c48bec549c (source withheld). Framingham Steatosis Index and albuminuria.
Spec: FSI per 1 unit (reading_1), albuminuria = ACR > 30 mg/g (reading_1), Model 3 (Table 2), NHANES 1999-2018,
fasting subsample weights (2/10) WTSAF4YR 1999-2002, (1/10) WTSAF2YR 2003-2018.
Sources in fulltext_txt/c48bec549c.txt: 'Data and study participants', 'Definition of FSI', 'Outcome ascertainment',
'Assessment of covariates', 'Statistical analysis', Table 1 (categories), Table 2 footnote.
"""
import numpy as np
import pandas as pd

REC_ID = "c48bec549c"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014",
          "2015-2016", "2017-2018"]


def _cat(df, name, conds, labels, base):
    df[name] = pd.Series(np.select(conds, labels, default=""), index=df.index).astype(object)
    df.loc[df[base].isna(), name] = np.nan
    return name


def _flag(pos, observed):
    return pd.Series(np.where(pos, 1.0, np.where(observed, 0.0, np.nan)), index=pos.index)


def _drink3(lib, df):
    """non / former / current drinker as defined in the paper. To 2016: nondrinker = < 12 drinks in lifetime
    (ALQ101/ALQ100/ALD100 = 2 and ALQ110 = 2); ever = >= 12 drinks in lifetime or any one year; current = ever and
    ALQ120Q > 0 (drank in past 12 months), former = ever and ALQ120Q = 0. 2017-2018 (no 12-drink items):
    nondrinker = ALQ111 2, former = ALQ111 1 and ALQ121 0, current = ALQ121 1-10."""
    S = lib.stems("ALQ")
    a101 = lib.raw(df, S, ["ALQ101", "ALQ100", "ALD100"])
    a110 = lib.raw(df, S, "ALQ110")
    q120 = lib.raw(df, S, "ALQ120Q").where(lambda v: v <= 365)
    ever = (a101 == 1) | (a110 == 1)
    old = np.select([(a101 == 2) & (a110 == 2), ever & (q120 > 0), ever & (q120 == 0)],
                    ["non", "current", "former"], default="")
    a111 = lib.raw(df, S, "ALQ111")
    a121 = lib.raw(df, S, "ALQ121")
    new = np.select([a111 == 2, a121.between(1, 10), (a111 == 1) & (a121 == 0)], ["non", "current", "former"], default="")
    isnew = df.cycle.isin(["2017-2018", "2017-2020"]).to_numpy()
    return pd.Series(np.where(isnew, new, old), index=df.index).replace("", np.nan).astype(object)


def _activity(lib, df):
    """vigorous yes/no and moderate yes/no. 1999-2006: PAD200 / PAD320 (vigorous / moderate leisure activity over past
    30 days; 1 yes, 2 no, 3 unable = no). 2007-2018 GPAQ: vigorous = PAQ605 or PAQ650 yes, moderate = PAQ620 or PAQ665
    yes; no when both items are no."""
    S = lib.stems("PAQ")
    vig, mod = pd.Series(np.nan, index=df.index), pd.Series(np.nan, index=df.index)
    p200, p320 = lib.raw(df, S, "PAD200"), lib.raw(df, S, "PAD320")
    old = df.cycle.isin(["1999-2000", "2001-2002", "2003-2004", "2005-2006"])
    vig[old] = p200[old].map({1: 1.0, 2: 0.0, 3: 0.0})
    mod[old] = p320[old].map({1: 1.0, 2: 0.0, 3: 0.0})
    a, b, c, d = (lib.raw(df, S, k) for k in ("PAQ605", "PAQ650", "PAQ620", "PAQ665"))
    v2 = _flag((a == 1) | (b == 1), (a == 2) & (b == 2))
    m2 = _flag((c == 1) | (d == 1), (c == 2) & (d == 2))
    vig[~old] = v2[~old]
    mod[~old] = m2[~old]
    return vig, mod


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["FSI"][0](df)
    df["y"] = (df.uacr > 30).astype(float).where(df.uacr.notna())
    S = lib.stems
    # covariate coding as Table 1 presents it
    _cat(df, "age2", [df.age < 60, df.age >= 60], ["lt60", "ge60"], "age")
    _cat(df, "bmi2", [df.bmi < 30, df.bmi >= 30], ["lt30", "ge30"], "bmi")
    _cat(df, "pir3", [df.pir < 1.3, df.pir < 3.5, df.pir >= 3.5], ["lt1.3", "1.3_3.5", "ge3.5"], "pir")
    _cat(df, "egfr2", [df.egfr_ckdepi2009 < 60, df.egfr_ckdepi2009 >= 60], ["lt60", "ge60"], "egfr_ckdepi2009")
    df["marital2"] = lib.raw(df, S("DEMO"), "DMDMARTL").map({1: "couple", 6: "couple", 2: "alone", 3: "alone",
                                                            4: "alone", 5: "alone"})
    # education for ages 18-19 from DMDEDUC3 (DMDEDUC2 covers 20+ only): grades up to 12 without diploma, 55, 66 = under
    # high school; 13 high school graduate, 14 GED = high school or equivalent; 15 more than high school
    e3 = lib.raw(df, S("DEMO"), "DMDEDUC3")
    e3c = pd.Series(np.where(e3.between(0, 12) | e3.isin([55, 66]), "lt_hs", np.where(e3.isin([13, 14]), "hs_ged",
                    np.where(e3 == 15, "gt_hs", ""))), index=df.index).replace("", np.nan)
    df["educ3"] = df.educ3.where(df.age >= 20, e3c)
    df["drink3"] = _drink3(lib, df)
    df["vig"], df["mod"] = _activity(lib, df)
    # diabetes: self-reported diagnosis, FPG >= 7.0 mmol/L (126 mg/dL), HbA1c >= 6.5 %, OGTT 2-h glucose >= 11.1 mmol/L (200)
    ogtt = lib.raw(df, S("OGTT"), "LBXGLT")
    dpos = (df.DIQ010 == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5) | (ogtt >= 200)
    df["dm"] = _flag(dpos, df.DIQ010.isin([1, 2, 3]) | df.glucose.notna() | df.hba1c.notna())
    # hypertension: self-reported, antihypertensive medication, SBP >= 130 or DBP >= 80 = hypertension_any_130
    # hyperlipidemia (NCEP): TG >= 150, TC >= 200, LDL-C >= 130, HDL-C <= 50 women / <= 40 men (mg/dL)
    lpos = ((df.tg >= 150) | (df.tc >= 200) | (df.ldl >= 130)
            | ((df.female == 1) & (df.hdl <= 50)) | ((df.female == 0) & (df.hdl <= 40)))
    df["hyperlip"] = _flag(lpos, df.tg.notna() & df.tc.notna() & df.hdl.notna())
    terms = ["x"]
    terms += lib.dummies(df, "age2", "lt60")
    terms += ["female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "bmi2", "lt30")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "marital2", "couple")
    terms += lib.dummies(df, "pir3", "lt1.3")
    terms += ["albumin_serum", "uric_acid", "hyperlip", "dm"]
    terms += lib.dummies(df, "drink3", "non")
    terms += ["hypertension_any_130", "vig", "mod"]
    terms += lib.dummies(df, "smoking", "never")
    terms += lib.dummies(df, "egfr2", "ge60")
    df["w"] = lib.pooled_weight(df, "fasting")
    dom = (df.age >= 18) & (df.pregnant != 1) & df.x.notna() & df.uacr.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "exposure scale of the Table 2 'Albuminuria' row": "reading_1: continuous FSI per 1 unit",
                "Albuminuria threshold": "reading_1: ACR > 30 mg/g",
                "FSI": "lib.INDEX['FSI'] (recon key; its hypertension and diabetes terms use the library composites "
                       "hypertension_any_140 and diabetes_any, not the covariate definitions of the paper)",
                "survey weight": "fasting weights pooled as stated ((2/10) WTSAF4YR, (1/10) WTSAF2YR)",
                "age, BMI, PIR, eGFR coding": "categories of Table 1 (methods give no coding): age < 60 / >= 60, "
                                              "BMI < 30 / >= 30, PIR < 1.3, 1.3-<3.5, >= 3.5, eGFR < 60 / >= 60",
                "eGFR equation": "not stated: CKD-EPI 2009 with race",
                "education": "educ3 (Table 1: under HS, HS or equivalent, college or above)",
                "marital": "married or living with partner (DMDMARTL 1, 6) vs living alone (2-5)",
                "race": "RIDRETH1 5 groups, ref non-Hispanic White",
                "smoking": "never/former/current (SMQ020, SMQ040) as defined",
                "alcohol": "non/former/current as defined; 2017-2018 approximated from ALQ111, ALQ121",
                "physical activity": "vigorous yes/no, moderate yes/no (PAD200/PAD320 1999-2006; GPAQ work or "
                                     "recreation items 2007-2018)",
                "diabetes": "DIQ010 1, FPG >= 126, HbA1c >= 6.5, OGTT 2-h >= 200 mg/dL (OGTT 2005-2016)",
                "hypertension": "hypertension_any_130 (BPQ020, BPQ050A, SBP >= 130, DBP >= 80)",
                "hyperlipidemia": "TG >= 150, TC >= 200, LDL >= 130, HDL <= 50 women / <= 40 men (mg/dL); no medication term",
                "albumin, uric acid": "continuous"},
            "notes": "domain = age >= 18, not pregnant, FSI computable, ACR present. Post-run technical fix 1: education "
                     "was missing for ages 18-19 (DMDEDUC2 is asked of adults 20+ only); added DMDEDUC3 for 18-19 "
                     "(wrong source column for that age band). Smoking and alcohol items of the adult SMQ/ALQ files "
                     "also miss part of ages 18-19; left as is (no youth file in the downloaded set)."}

"""M1 builder, 988df86a79 (source withheld). DII per SD, CKD, Model 3.
Spec: NHANES 2011-2018, age >= 18, CKD = eGFR < 60 (CKD-EPI, version not stated) or uACR >= 30 mg/g, no weighting
reported (weight_variable null), Table 2 Model 3 OR 1.16 per SD.
Sources (fulltext_txt/988df86a79.txt): exclusions lines 63-66, data collection and definitions lines 68-104,
DII lines 107-121, models lines 128-134, Table 1 coding lines 215-370.
"""
import numpy as np, pandas as pd

REC_ID = "988df86a79"
CYCLES = ["2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def _clean(s, bad):
    return s.where(~s.isin(bad))


def _drinks_per_day(lib, df):
    """ALQ130 average drinks per day on drinking days in the past 12 months; 0 for participants skipped because they
    did not drink in the past year (ALQ120Q = 0 or ALQ121 = 0) or are lifetime non-drinkers (ALQ110 = 2 / ALQ111 = 2)."""
    s = lib.stems("ALQ")
    a130 = _clean(lib.raw(df, s, "ALQ130"), [777, 999])
    q = _clean(lib.raw(df, s, "ALQ120Q"), [777, 999])
    a121 = _clean(lib.raw(df, s, "ALQ121"), [77, 99])
    yr12 = _clean(lib.raw(df, s, "ALQ101"), [7, 9])
    life12 = _clean(lib.raw(df, s, "ALQ110"), [7, 9])
    a111 = _clean(lib.raw(df, s, "ALQ111"), [7, 9])
    zero = (q == 0) | (a121 == 0) | ((yr12 == 2) & (life12 == 2)) | (a111 == 2)
    return a130.where(~(zero & a130.isna()), 0.0)


def build(lib):
    df = lib.table(CYCLES)
    df["dii"] = lib.DIET[REC_ID](df, reading=1)
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())
    dom = (df.age >= 18) & df.dii.notna() & df.uacr.notna() & e.notna()
    df["x"] = lib.sd_scale(df, "dii", dom)
    # current smoking (yes vs no): SMQ040 every day / some days among SMQ020 = 1; never smokers 0
    df["cur_smoke"] = pd.Series(np.where((df.SMQ020 == 1) & df.SMQ040.isin([1, 2]), 1.0,
                                np.where((df.SMQ020 == 2) | ((df.SMQ020 == 1) & (df.SMQ040 == 3)), 0.0, np.nan)), index=df.index)
    df["alc_dpd"] = _drinks_per_day(lib, df)
    df["albumin_gL"] = df.albumin_serum * 10
    # medication covariates from the questionnaires (BPQ050A, BPQ100D, DIQ050/DIQ070); not asked (skip) = 0
    bq = lib.stems("BPQ")
    bp_rec = df.BPQ020.isin([1, 2])
    a = lib.raw(df, bq, "BPQ050A")
    df["antihyp"] = pd.Series(np.where(a == 1, 1.0, np.where(bp_rec & ~a.isin([7, 9]), 0.0, np.nan)), index=df.index)
    l = lib.raw(df, bq, "BPQ100D")
    df["lipid_drug"] = pd.Series(np.where(l == 1, 1.0, np.where(bp_rec & ~l.isin([7, 9]), 0.0, np.nan)), index=df.index)
    di_rec = df.DIQ010.isin([1, 2, 3])
    df["glu_drug"] = pd.Series(np.where((df.DIQ050 == 1) | (df.DIQ070 == 1), 1.0, np.where(di_rec, 0.0, np.nan)), index=df.index)
    # diabetes: self-reported diagnosis or FPG >= 7.0 mmol/L (126 mg/dL) or glucose-lowering drugs
    pos = (df.DIQ010 == 1) | (df.glucose >= 126) | (df.glu_drug == 1)
    df["dm"] = pd.Series(np.where(pos, 1.0, np.where(di_rec | df.glucose.notna(), 0.0, np.nan)), index=df.index)
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["pir", "bmi", "sbp", "dbp", "cur_smoke", "alc_dpd", "glucose", "albumin_gL", "tg", "tc", "hdl", "ldl",
              "dm", "antihyp", "lipid_drug", "glu_drug"]
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": None,
            "readings": {
                "n_analytic (spec ambiguity)": "reading_1: 19,317 (no effect on the builder)",
                "DII": "diet_indices_i139 DIET['988df86a79'] reading 1: 28 named parameters, day-1 recall (stated)",
                "per SD": "DII / SD of DII over the domain rows (unweighted)",
                "weighting": "none reported (Empower, 'multivariable logistic regression'); unweighted, weight None",
                "age": "'Age under 18' excluded read as age >= 18 kept",
                "dietary questionnaire and blood samples": "expressed by non-missing DII, uACR and eGFR",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated)",
                "CKD threshold": "uACR >= 30 mg/g (reading_1 of the spec threshold)",
                "race": "RIDRETH1 5 groups as in Table 1 (Other races includes Asian)",
                "current smoking": "binary current smoker (SMQ040 1/2) vs not (never or former), as in Table 1 'Current smoking'",
                "alcohol intake": "drinks per day ALQ130 (Table 1 'drinks per day'); 0 for past-year and lifetime non-drinkers who were skipped",
                "FPG": "fasting plasma glucose LBXGLU (fasting subsample); TG and LDL-C from the fasting TRIGLY file",
                "albumin": "serum albumin g/L (LBXSAL x 10)",
                "diabetes mellitus": "DIQ010 = 1 or FPG >= 126 mg/dL or DIQ050/DIQ070 = 1 (borderline = no)",
                "antihypertensive drugs": "BPQ050A = 1 taking prescribed medicine for high blood pressure; skipped = 0",
                "lipoprotein-lowering drugs": "BPQ100D = 1 taking prescribed medicine for cholesterol; skipped = 0",
                "glucose-lowering drugs": "DIQ050 insulin or DIQ070 diabetic pills = 1; skipped = 0",
                "PIR, BMI, SBP, DBP, TC, HDL": "continuous",
            },
            "notes": "written before first execution; no post-run changes"}

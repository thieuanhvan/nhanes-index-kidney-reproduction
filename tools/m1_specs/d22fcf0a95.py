"""M1 builder, d22fcf0a95 (NPAR and CKD, PLoS One 2024, doi [DOI withheld]).
Spec: NPAR quartiles, Model 3, Q2 vs Q1 (reading_1 of both ambiguities), CKD = eGFR < 60 or UACR >= 30, NHANES 2009-2018.
Sources (fulltext_txt/d22fcf0a95.txt): population line 12, outcome line 20, covariate coding line 24,
models and weighting line 27, Table 1 lines 104-165, Table 2 lines 168-173.
"""
import numpy as np
import pandas as pd

REC_ID = "d22fcf0a95"
CYCLES = ["2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def labels(index, pairs):
    """object Series; first matching (mask, label) wins; NaN when no mask is True."""
    out = pd.Series(np.nan, index=index, dtype=object)
    for m, lab in reversed(pairs):
        out[m.fillna(False).to_numpy(bool)] = lab
    return out


ALQ121_DAYS = {0: 0.0, 1: 365.0, 2: 286.0, 3: 182.0, 4: 104.0, 5: 52.0, 6: 30.0, 7: 12.0, 8: 9.0, 9: 4.5, 10: 1.5}


def alcohol3(df, lib):
    """never (<12 drinks in lifetime) / low-to-moderate / severe (> 1 drink/day women, > 2 men, past-year average)."""
    st = lib.stems("ALQ")
    a110 = lib.raw(df, st, "ALQ110")
    a111 = lib.raw(df, st, "ALQ111")
    q = lib.raw(df, st, "ALQ120Q").where(lambda v: v <= 365)
    u = lib.raw(df, st, "ALQ120U").map({1: 52.0, 2: 12.0, 3: 1.0})
    f121 = lib.raw(df, st, "ALQ121").map(ALQ121_DAYS)
    d130 = lib.raw(df, st, "ALQ130").where(lambda v: v < 77)
    days = (q * u).where(q != 0, 0.0).where(q.notna())
    days = days.fillna(f121)
    avg = (days / 365.0 * d130).where(days != 0, 0.0).where(days.notna())
    never = (a110 == 2) | (a111 == 2)
    lim = np.where(df.female == 1, 1.0, 2.0)
    return labels(df.index, [(never, "never"), (avg > lim, "severe"), (avg.notna(), "lowmod")])


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["NPAR"][0](df)
    # CKD: eGFR (CKD-EPI 2009 with race, version not stated) < 60 or UACR >= 30 mg/g
    e = df.egfr_ckdepi2009
    df["y"] = ((e < 60) | (df.uacr >= 30)).astype(float).where(e.notna() & df.uacr.notna())
    dom = (df.age >= 20) & (df.pregnant != 1) & df.scr.notna() & df.uacr.notna() & df.x.notna()
    # exposure: NPAR quartiles cut on the analytic population (unweighted; Table 1 quartile n are equal)
    qn, cuts = lib.quantile_dummies(df, "x", dom, k=4)
    # covariates (line 24 coding)
    mar = lib.raw(df, lib.stems("DEMO"), "DMDMARTL")
    df["marital6"] = mar.map({1: "married", 2: "widowed", 3: "divorced", 4: "separated", 5: "never", 6: "partner"})
    df["pir3"] = labels(df.index, [(df.pir < 1.5, "low"), (df.pir <= 3.5, "mid"), (df.pir > 3.5, "high")])
    df["alc3"] = alcohol3(df, lib)
    bpmed = lib.raw(df, lib.stems("BPQ"), "BPQ050A")
    df["hbp"] = flag((df.sbp >= 130) | (df.dbp >= 90) | (df.dx_hypertension == 1) | (bpmed == 1),
                     df.dx_hypertension.notna() | (df.sbp.notna() & df.dbp.notna()))
    pills = lib.raw(df, lib.stems("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, lib.stems("DIQ"), "DIQ050")
    ogtt = lib.raw(df, lib.stems("OGTT"), "LBXGLT")
    df["dm"] = flag((df.dx_diabetes == 1) | (pills == 1) | (ins == 1) | (df.glucose >= 126) | (ogtt >= 200) | (df.hba1c >= 6.5),
                    df.dx_diabetes.notna() | df.hba1c.notna() | df.glucose.notna())
    mcq = {k: lib.raw(df, lib.stems("MCQ"), k) for k in ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E"]}
    m = pd.DataFrame(mcq)
    df["cvd4"] = flag((m == 1).any(axis=1), m.isin([1, 2]).all(axis=1))
    terms = list(qn)
    terms += ["age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "marital6", "married")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "pir3", "low")
    terms += lib.dummies(df, "alc3", "never")
    terms += lib.dummies(df, "smoking", "never")
    terms += ["bmi", "hbp", "dm", "cvd4"]
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x_q2", "weight": None,
            "readings": {
                "estimate (spec ambiguity 1 and 2)": "reading_1: Model 3, NPAR Q2 vs Q1",
                "survey weighting": "unweighted (no weight named; 'according to the NHANES recommendations' only; Table 1 reports raw n and unweighted %)",
                "quartile cut": "unweighted quartiles of NPAR in the analytic population (age >= 20, not pregnant, Scr, UACR and NPAR present)",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated)",
                "UACR threshold": ">= 30 mg/g ('30 mg/g or higher')",
                "age, BMI": "continuous (Table 1 mean +- SD, no categories stated)",
                "race": "RIDRETH1 5 groups, reference non-Hispanic white",
                "marital status": "DMDMARTL 6 categories as listed, reference married",
                "education": "3 categories (DMDEDUC2 1-2 / 3 / 4-5)",
                "PIR": "low < 1.5, medium 1.5 to 3.5, high > 3.5",
                "alcohol": "never = ALQ110 = 2 (< 12 drinks in lifetime; 2017-2018 ALQ111 = 2 never had a drink, the only lifetime item); past-year average drinks/day = days drunk per year (ALQ120Q x ALQ120U, or ALQ121 class midpoints in 2017-2018) / 365 x ALQ130; severe > 1 (women) or > 2 (men); others low-to-moderate",
                "smoking": "harmonised never/former/current (SMQ020, SMQ040)",
                "HBP": "SBP >= 130 or DBP >= 90 (as written), BPQ020 = 1, or BPQ050A = 1",
                "DM": "DIQ010 = 1, pills DIQ070 = 1, insulin DIQ050, fasting glucose >= 126 mg/dL (7.0 mmol/L), OGTT 2 h LBXGLT >= 200 mg/dL (11.1 mmol/L; no OGTT file in 2017-2018), HbA1c >= 6.5",
                "CVD": "MCQ160B, C, D, E (heart failure, CHD, angina, heart attack; stroke not listed)",
                "missing covariates": "paper imputed (random forest MICE); here complete case on model terms as the M1 fit does"},
            "notes": f"NPAR quartile cut points on domain: {[float(np.trunc(c * 1e4) / 1e4) for c in cuts]} (truncated)"}

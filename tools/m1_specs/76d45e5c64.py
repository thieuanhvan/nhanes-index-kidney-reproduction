"""M1 builder, 76d45e5c64 (UHR and DKD; source withheld).
Spec: UHR (uric acid / HDL-C x 100) Q4 vs Q1, DKD = UACR > 30 mg/g and/or eGFR (CKD-EPI 2009) < 60 in diabetes,
Model 3, NHANES 2001-2018. Sources in fulltext_txt/76d45e5c64.txt: population line 12, exposure and outcome line 16,
covariates lines 19-20, Model 3 line 23, covariate denominators Table 1 (line 41).
"""
import numpy as np, pandas as pd
import diet_indices_i139 as di

REC_ID = "76d45e5c64"
CYCLES = ["2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016",
          "2017-2018"]
ALQ142_DPM = {1: 30.0, 2: 25.0, 3: 15.0, 4: 8.6, 5: 4.3, 6: 2.5, 7: 1.0, 8: 0.75, 9: 0.375, 10: 0.125}   # days per month


def _flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def _bio(lib, df, col):
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for c in df.cycle.unique():
        x = lib.nc.read(c, "BIOPRO")
        if x is None or col not in x.columns:
            continue
        m = (df.cycle == c).to_numpy()
        out[m] = x.drop_duplicates("SEQN").set_index("SEQN")[col].reindex(df.loc[m, "SEQN"]).to_numpy(float)
    return out


def _per_month(q, u):
    q = q.where(q.between(0, 366))
    return q * u.map({1.0: 30 / 7, 2.0: 1.0, 3.0: 1 / 12})


def _drinking(lib, df):
    """mild / moderate / heavy among current drinkers (past 12 months); never and former drinkers missing.
    heavy: >= 3 drinks/day women, >= 4 men, or binge >= 5 days/month; moderate: >= 2 women, >= 3 men, or binge >= 2 days/month."""
    s = lib.stems("ALQ")
    freq = lib.raw(df, s, "ALQ120Q"); f17 = lib.raw(df, s, "ALQ121")
    j = (df.cycle == "2017-2018")
    current = pd.Series(np.where(j, np.where(f17.between(1, 10), 1.0, np.where(f17 == 0, 0.0, np.nan)),
                                 np.where(freq.between(1, 366), 1.0, np.where(freq == 0, 0.0, np.nan))), index=df.index)
    dpd = lib.raw(df, s, "ALQ130").where(lambda v: v.between(1, 76))
    b140 = _per_month(lib.raw(df, s, "ALQ140Q"), lib.raw(df, s, "ALQ140U"))
    b141 = _per_month(lib.raw(df, s, "ALQ141Q"), lib.raw(df, s, "ALQ141U"))
    b142 = lib.raw(df, s, "ALQ142").map(lambda v: 0.0 if v == 0 else ALQ142_DPM.get(v, np.nan))
    binge = b140.fillna(b141).fillna(b142)
    f = df.female == 1
    heavy = (f & (dpd >= 3)) | (~f & (dpd >= 4)) | (binge >= 5)
    moder = (f & (dpd >= 2)) | (~f & (dpd >= 3)) | (binge >= 2)
    cat = pd.Series(np.where(heavy, "heavy", np.where(moder, "moderate", "mild")), index=df.index).astype(object)
    return cat.where((current == 1) & dpd.notna(), np.nan)


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["UHR_x100"][0](df)
    egfr = df.egfr_ckdepi2009
    df["y"] = ((df.uacr > 30) | (egfr < 60)).astype(float).where(df.uacr.notna() & egfr.notna())

    pills = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    glu_si = _bio(lib, df, "LBDSGLSI")
    pos = ((df.DIQ010 == 1) | (df.DIQ050 == 1) | (pills == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5)
           | (glu_si >= 11.1))
    obs = df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna() | glu_si.notna()
    df["dm"] = _flag(pos, obs)

    mar = lib.raw(df, s("DEMO"), "DMDMARTL")
    df["coupled"] = mar.map({1.0: 1.0, 6.0: 1.0, 2.0: 0.0, 3.0: 0.0, 4.0: 0.0, 5.0: 0.0})
    df["drink"] = _drinking(lib, df)
    pa = di.lifestyle(df).pa_total
    df["pa_high"] = (pa >= 500).astype(float).where(pa.notna())
    ins = df.DIQ050.map({1.0: 1.0, 2.0: 0.0}); pil = pills.map({1.0: 1.0, 2.0: 0.0})
    df["antidm"] = _flag((ins == 1) | (pil == 1), ins.notna() | pil.notna())
    df["antilip"] = lib.raw(df, s("BPQ"), "BPQ100D").map({1.0: 1.0, 2.0: 0.0})
    g = lib.raw(df, s("DIQ"), "DID040G"); q = lib.raw(df, s("DIQ"), "DID040Q"); a = lib.raw(df, s("DIQ"), "DID040")
    age_dx = pd.Series(np.where(g == 1, q, np.where(g == 2, 0.0, a)), index=df.index)
    age_dx = age_dx.where(age_dx != 666, 0.0)
    age_dx = age_dx.where(age_dx < 100)
    df["dm_dur"] = (df.age - age_dx).clip(lower=0).where(df.DIQ010 == 1)
    df["htn"] = _flag((df.BPQ020 == 1) | (df.BPQ_med_now == 1) | (df.sbp >= 130) | (df.dbp >= 80),
                      df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna()))
    df["bmi_g"] = np.select([df.bmi < 25, df.bmi < 30, df.bmi >= 30], ["lt25", "25_30", "ge30"], default="").astype(object)
    df.loc[df.bmi.isna(), "bmi_g"] = np.nan

    dom = (df.age >= 20) & (df.dm == 1) & df.x.notna() & df.uacr.notna() & egfr.notna()
    qn, cuts = lib.quantile_dummies(df, "x", dom, 4, prefix="x_q")
    terms = qn + ["age", "female"]
    terms += lib.dummies(df, "race1", "mexican_american")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["pir", "coupled"]
    terms += lib.dummies(df, "smoking", "never")
    terms += lib.dummies(df, "drink", "mild")
    terms += ["pa_high", "antidm", "antilip", "dm_dur", "htn", "cvd_history"]
    terms += lib.dummies(df, "bmi_g", "lt25")
    terms += ["glucose", "hba1c", "ldl", "wbc", "hgb", "neut_n"]

    df["w"] = lib.pooled_weight(df, "fasting")
    readings = {
        "exposure scale": "Q4 vs Q1 of UHR (reading_1); quartiles cut unweighted on the domain (rule 6); the paper's printed cuts are 9.27, 12.27, 16.11",
        "eGFR": "CKD-EPI 2009 with race (stated)",
        "UACR threshold": "> 30 mg/g as written",
        "diabetes": "DIQ010=1, insulin DIQ050=1, pills (DIQ070, DID070 2005-2008), FPG >= 126 mg/dL, HbA1c >= 6.5, serum glucose LBDSGLSI >= 11.1 mmol/L as random glucose; borderline not diabetes",
        "pregnancy": "not excluded (not listed by the paper)",
        "survey weighting": "weighted, fasting subsample weight (Model 3 includes fasting plasma glucose and LDL-C, both fasting subsample)",
        "drinking": "mild/moderate/heavy among past-12-month drinkers by the usual NHANES rule (drinks per drinking day ALQ130 and binge days per month from ALQ140/ALQ141/ALQ142); never and former drinkers have no category in the paper and are missing (Table 1 denominators 3,479 of 7,138)",
        "physical activity": "< 500 vs >= 500 MET-min/week (Table 1), MET-min/week from tools/diet_indices_i139.lifestyle (GPAQ 2007-2018; leisure PAQIAF + transport 2001-2006)",
        "anti-diabetic agents": "insulin DIQ050 or pills DIQ070/DID070 yes; no when answered no; not asked = missing (Table 1 denominator 5,635)",
        "anti-hyperlipidemia agents": "BPQ100D now taking prescribed medicine for cholesterol (1 yes, 2 no); not asked = missing (Table 1 denominator 3,536)",
        "diabetic duration": "age minus age at diagnosis (DID040G/DID040Q 2001-2004, DID040 2005+), diagnosed diabetes only (Table 1 reports a duration only for those asked)",
        "hypertension": "BPQ020=1 or BPQ050A=1 or mean SBP >= 130 or DBP >= 80 (2017 ACC/AHA)",
        "CVD": "self-reported heart failure, CHD, angina, MI or stroke (cvd_history)",
        "covariate coding": "race 5 groups, education 3 levels, marital coupled (married, living with partner) vs single, smoking never/former/current, BMI <25/25-30/>=30, age, PIR, FPG, HbA1c, LDL, WBC, hemoglobin, neutrophils continuous",
    }
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x_q4", "weight": "w",
            "readings": readings, "notes": "builder written by an M1 agent 15/11/2026; quartile cuts on the domain " +
            ", ".join(str(int(c * 100) / 100) for c in cuts)}

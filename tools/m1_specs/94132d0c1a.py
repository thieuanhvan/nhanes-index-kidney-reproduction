"""M1 builder, 94132d0c1a (SIRI and CKD, Front Endocrinol 2024, doi [DOI withheld]).
Spec: SIRI per 1 unit, CKD = ACR >= 30 mg/g or eGFR (CKD-EPI 2009) < 60, Model 3, NHANES 1999-2016 + 2017-March 2020.
Sources in fulltext_txt/94132d0c1a.txt: population line 16, CKD line 20, covariates and HTN/DM definitions lines 23-24,
Model 3 line 27 and Table 2 footnote, covariate coding Table 1 (lines 95-157).
"""
import numpy as np, pandas as pd

REC_ID = "94132d0c1a"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014",
          "2015-2016", "2017-2020"]


def _bio(lib, df, cols):
    """Column of the standard biochemistry file (LAB18 / L40_B / L40_C / BIOPRO_*) for every cycle."""
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for c in df.cycle.unique():
        x = lib.nc.read(c, "BIOPRO")
        if x is None:
            continue
        col = next((k for k in cols if k in x.columns), None)
        if col is None:
            continue
        m = (df.cycle == c).to_numpy()
        out[m] = x.drop_duplicates("SEQN").set_index("SEQN")[col].reindex(df.loc[m, "SEQN"]).to_numpy(float)
    return out


def _yn(v):
    """1 yes / 2 no questionnaire item to 1/0, other codes missing."""
    return v.map({1.0: 1.0, 2.0: 0.0})


def _cat(df, name, conds, labels):
    df[name] = np.select(conds, labels, default="").astype(object)
    df.loc[df[name] == "", name] = np.nan


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["SIRI"][0](df)
    egfr = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (egfr < 60)).astype(float).where(df.uacr.notna() & egfr.notna())

    # covariates coded as Table 1 presents them
    _cat(df, "age_g", [df.age <= 40, df.age <= 60, df.age > 60], ["20_40", "41_60", "gt60"])
    df.loc[df.age.isna(), "age_g"] = np.nan
    _cat(df, "bmi_g", [df.bmi < 25, df.bmi < 30, df.bmi >= 30], ["normal", "over", "obese"])
    df.loc[df.bmi.isna(), "bmi_g"] = np.nan
    _cat(df, "pir_g", [df.pir < 1.3, df.pir < 3.5, df.pir >= 3.5], ["low", "medium", "high"])
    df.loc[df.pir.isna(), "pir_g"] = np.nan
    df["smk100"] = _yn(df.SMQ020)
    alc = lib.raw(df, s("ALQ"), ["ALQ100", "ALD100", "ALQ101", "ALQ111"])
    df["alcohol"] = _yn(alc)
    mar = lib.raw(df, s("DEMO"), ["DMDMARTL", "DMDMARTZ"])
    p = df.cycle == "2017-2020"
    m_l = mar.map({1.0: "married", 5.0: "never", 6.0: "partner", 2.0: "other", 3.0: "other", 4.0: "other"})
    m_z = mar.map({1.0: "married", 2.0: "other", 3.0: "never"})
    df["marital"] = m_l.where(~p, m_z).astype(object)

    # diabetes (line 24): self-report, insulin or prescription drugs, HbA1c > 6.5, fasting glucose >= 7.0 mmol/L
    pills = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    pos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (pills == 1) | (df.hba1c > 6.5) | (df.glucose >= 126)
    obs = df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna()
    df["dm"] = np.where(pos, 1.0, np.where(obs, 0.0, np.nan))
    # hypertension (line 24): self-report, mean BP above 130/80 ("surpassed"), taking prescription
    pos = (df.BPQ020 == 1) | (df.BPQ_med_now == 1) | (df.sbp > 130) | (df.dbp > 80)
    obs = df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna())
    df["htn"] = np.where(pos, 1.0, np.where(obs, 0.0, np.nan))

    df["phos"] = _bio(lib, df, ["LBXSPH", "LBDSPH"])
    df["tg_s"] = df.tg_serum

    terms = ["x", "female"]
    terms += lib.dummies(df, "age_g", "20_40")
    terms += lib.dummies(df, "race1", "mexican_american")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "bmi_g", "normal")
    terms += ["smk100", "alcohol", "uric_acid", "phos", "tc", "tg_s", "ast", "alt"]
    terms += lib.dummies(df, "pir_g", "low")
    terms += ["cvd_history"]
    terms += lib.dummies(df, "marital", "married")
    terms += ["dm", "htn"]

    df["w"] = lib.pooled_weight(df, "mec")
    cancer = lib.raw(df, s("MCQ"), "MCQ220") == 1
    dom = (df.age >= 20) & (df.pregnant != 1) & ~cancer & df.uacr.notna() & egfr.notna() & df.x.notna()
    readings = {
        "UACR threshold": "ACR >= 30 mg/g (reading_1)",
        "cycles": "1999-2000 to 2015-2016 plus the 2017-March 2020 pre-pandemic files, pooled with NCHS weights (reading_1)",
        "survey weighting": "weighted, MEC weight (CBC, biochemistry, urine all MEC components; no fasting restriction in the paper's flow)",
        "eGFR": "CKD-EPI 2009 with race (stated)",
        "triglycerides": "serum biochemistry TG (LBXSTR): the paper's flow (line 16) lists no fasting restriction, so its 41,089 include non-fasting participants",
        "covariate coding": "as Table 1 presents: age 20-40/41-60/>60, BMI <25/25-29.9/>=30, PIR low/medium/high with cuts 1.3 and 3.5 (cut points not given; common NHANES cuts), smoking >=100 cigarettes yes/no, alcohol yes/no, marital married/never/partner/other, education 3 levels",
        "alcohol": "yes = had at least 12 drinks in any one year (ALQ100/ALD100/ALQ101) 1999-2016; ever had a drink (ALQ111) in 2017-2020 (the 12-drink item is not asked there)",
        "marital 2017-2020": "DMDMARTZ 1 (married or living with partner) coded married, 2 other, 3 never married",
        "education Others": "Table 1 'Others' (n 42) are refused/don't know codes, kept missing (rule 7)",
        "diabetes": "self-report DIQ010=1 or insulin DIQ050=1 or pills (DIQ070, DID070 in 2005-2008) or HbA1c > 6.5 (strict, as written) or fasting glucose >= 126 mg/dL; borderline not diabetes",
        "hypertension": "BPQ020=1 or BPQ050A=1 or mean SBP > 130 or mean DBP > 80 ('surpassed')",
        "CVD": "self-reported CHF, CHD, angina, MI or stroke (not defined in the paper; cvd_history)",
        "cancer": "MCQ220 = 1 excluded; missing MCQ220 kept",
        "age, sex, race": "race RIDRETH1 5 groups as Table 1",
    }
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": readings, "notes": "builder written by an M1 agent 15/11/2026"}

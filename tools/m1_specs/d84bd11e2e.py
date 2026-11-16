"""M1 builder, d84bd11e2e (J Health Popul Nutr 2025, doi [DOI withheld]). CMI and DKD in diabetes.
Spec: Log2-CMI per 1 log2 unit, DKD = eGFR < 60 and/or UACR >= 30 mg/g among participants with diabetes (reading_1),
Model 3 (Table 3; reading_1: no BMI), NHANES 1999-2018, weighted.
Sources in fulltext_txt/d84bd11e2e.txt: 'Study population and participant selection', Table 1 (CMI formula, mmol/L),
'Outcome variables', 'Covariates', 'Statistical analysis', Table 3.
"""
import numpy as np
import pandas as pd

REC_ID = "d84bd11e2e"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014",
          "2015-2016", "2017-2018"]


def _cat(df, name, conds, labels, base):
    df[name] = pd.Series(np.select(conds, labels, default=""), index=df.index).astype(object)
    df.loc[df[base].isna(), name] = np.nan
    return name


def _flag(pos, observed):
    return pd.Series(np.where(pos, 1.0, np.where(observed, 0.0, np.nan)), index=pos.index)


def build(lib):
    from diet_indices_i139 import _paq   # read-only use of the shared MET-min/week helper (GPAQ 2007-2018, PAQIAF 1999-2006)
    df = lib.table(CYCLES)
    cmi = lib.INDEX["CMI"][0](df)
    df["x"] = np.log2(cmi.where(cmi > 0))
    ok = df.uacr.notna() & df.egfr_ckdepi2009.notna()
    df["y"] = ((df.uacr >= 30) | (df.egfr_ckdepi2009 < 60)).astype(float).where(ok)
    S = lib.stems
    # diabetes: self-report, insulin or other diabetes medication, FPG >= 7.0 mmol/L (126 mg/dL), HbA1c >= 6.5
    pills = lib.raw(df, S("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, S("DIQ"), "DIQ050")
    dpos = (df.DIQ010 == 1) | (ins == 1) | (pills == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5)
    df["dm"] = _flag(dpos, df.DIQ010.isin([1, 2, 3]) | df.glucose.notna() | df.hba1c.notna())
    df["race4"] = df.race1.map({"nh_white": "nh_white", "nh_black": "nh_black", "mexican_american": "mexican_american",
                                "other_hispanic": "other", "other_multi": "other"})
    _cat(df, "pir3", [df.pir <= 1.39, df.pir <= 3.49, df.pir > 3.49], ["low", "medium", "high"], "pir")
    pa = _paq(df)
    df["pa_met"] = pd.to_numeric(pa["pa_total"], errors="coerce")
    df["pa_active"] = (df.pa_met >= 500).astype(float).where(df.pa_met.notna())
    # hyperlipidemia: HDL <= 40 men / <= 50 women, TG >= 150, TC >= 200, LDL >= 130 mg/dL, or cholesterol-lowering drug
    b100d = lib.raw(df, S("BPQ"), "BPQ100D")
    lpos = ((df.tg >= 150) | (df.tc >= 200) | (df.ldl >= 130) | (b100d == 1)
            | ((df.female == 1) & (df.hdl <= 50)) | ((df.female == 0) & (df.hdl <= 40)))
    df["hyperlip"] = _flag(lpos, df.tg.notna() & df.tc.notna() & df.hdl.notna())
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race4", "nh_white")
    terms += lib.dummies(df, "pir3", "low")
    terms += ["pa_active"]
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "smoking", "never")
    terms += ["cvd_history", "hypertension_any_140", "hyperlip"]
    df["w"] = lib.pooled_weight(df, "fasting")
    dom = ((df.age >= 20) & (df.pregnant != 1) & (df.dm == 1) & df.x.notna() & df.bmi.notna() & df.y.notna())
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "BMI in Model 2/Model 3 adjustment set": "reading_1: Table 3 footnote set, no BMI",
                "population / outcome": "reading_1: DKD (eGFR < 60 and/or UACR >= 30) among participants with diabetes "
                                        "by the paper's criteria, non-diabetes excluded",
                "diabetes": "DIQ010 1, insulin (DIQ050), pills (DIQ070; DID070 2005-2008), FPG >= 126, HbA1c >= 6.5",
                "eGFR equation": "not stated: CKD-EPI 2009 with race",
                "exposure": "log2 of lib.INDEX['CMI'] (mmol/L as stated); a unit constant only shifts log2-CMI",
                "index exclusions": "missing data for CMI, ABSI, BRI, VAI, LAP: CMI computable and BMI present (waist, "
                                    "height, TG, HDL-C, BMI)",
                "survey weight": "fasting subsample weight (CMI needs fasting TG; variable not named)",
                "age": "continuous years",
                "race": "4 groups as listed: NH White (ref), NH Black, Mexican American, other (other Hispanic + other)",
                "PIR": "methods cut points: <= 1.39, > 1.39 to <= 3.49, > 3.49",
                "physical activity": "MET-min/week >= 500 active (diet_indices_i139._paq 'pa_total': GPAQ 2007-2018 "
                                     "work + transport + recreation; 1999-2006 leisure PAQIAF + transport, MET 8/4)",
                "education": "educ3 (less than HS, HS graduate/GED, above HS) as the methods list it",
                "smoking": "never/former/current",
                "CVD": "cvd_history: CHF, CHD, angina, heart attack, stroke (MCQ160B-F)",
                "hypertension": "hypertension_any_140 (BPQ020, BPQ050A, SBP >= 140, DBP >= 90)",
                "hyperlipidemia": "HDL <= 40 men / <= 50 women, TG >= 150, TC >= 200, LDL >= 130 mg/dL, BPQ100D 1"},
            "notes": "marital status is listed among covariates but is in none of the model sets; not adjusted."}

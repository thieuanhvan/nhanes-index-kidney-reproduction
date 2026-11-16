"""M1 builder, edaf91172c (HRR and CKD, Medicine 2024, doi [DOI withheld]).
Spec: HRR continuous per 1 unit, Model 3 (Table 2 footnote, 14 covariates), CKD = eGFR < 60 or UACR >= 30, NHANES 2005-2018, weighted.
Sources (fulltext_txt/edaf91172c.txt): population (Methods 2.1), outcome (2.3), covariates (2.4), weights and models (2.5),
Table 1 coding, Table 2.
"""
import numpy as np
import pandas as pd

REC_ID = "edaf91172c"
CYCLES = ["2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def labels(index, pairs):
    out = pd.Series(np.nan, index=index, dtype=object)
    for m, lab in reversed(pairs):
        out[m.fillna(False).to_numpy(bool)] = lab
    return out


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["HRR"][0](df)
    e = df.egfr_ckdepi2009
    df["y"] = ((e < 60) | (df.uacr >= 30)).astype(float).where(e.notna() & df.uacr.notna())
    dom = (df.age >= 20) & df.scr.notna() & df.uacr.notna() & df.x.notna()
    df["age50"] = flag(df.age >= 50, df.age.notna())
    df["race4"] = df.race1.map({"mexican_american": "mex", "nh_white": "white", "nh_black": "black",
                                "other_hispanic": "other", "other_multi": "other"})
    df["bmi3"] = labels(df.index, [(df.bmi < 25, "normal"), (df.bmi < 30, "over"), (df.bmi >= 30, "obese")])
    df["pir3"] = labels(df.index, [(df.pir < 1, "lt1"), (df.pir <= 3, "1to3"), (df.pir > 3, "gt3")])
    smq = lib.raw(df, lib.stems("SMQ"), "SMQ020")
    df["smoke"] = flag(smq == 1, smq.isin([1, 2]))
    df["htn"] = df["hypertension_any_130"]
    hc = lib.raw(df, lib.stems("BPQ"), "BPQ080")
    df["hichol"] = flag(hc == 1, hc.isin([1, 2]))
    pills = lib.raw(df, lib.stems("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, lib.stems("DIQ"), "DIQ050")
    df["dm"] = flag((df.dx_diabetes == 1) | (df.hba1c > 6.5) | (df.glucose_serum >= 200) | (pills == 1) | (ins == 1),
                    df.dx_diabetes.notna() | df.hba1c.notna() | df.glucose_serum.notna())
    terms = ["x", "age50", "female"]
    terms += lib.dummies(df, "bmi3", "normal")
    terms += lib.dummies(df, "race4", "white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "pir3", "lt1")
    terms += ["smoke", "htn", "hichol", "cvd_history", "dm", "bun", "alt", "ast"]
    df["w"] = lib.pooled_weight(df, "mec")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "scale (spec ambiguity 1)": "reading_1: OR per 1-unit HRR",
                "covariate list (spec ambiguity 2)": "reading_1: the 14 variables of the Table 2 footnote",
                "survey weighting": "weighted ('weights per official NHANES documentation'), MEC weight (CBC and biochemistry), 2-year weight x 2 / 14",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated)",
                "UACR threshold": ">= 30 mg/g",
                "age": "as coded in Table 1: < 50 vs >= 50 (Table 1 shows '<50' and '>50'; age 50 put in the upper group); methods give no coding",
                "BMI": "as coded in Table 1: normal < 25 (reference), overweight 25 to < 30, obese >= 30",
                "race": "4 groups as in Table 1: Mexican American, NH black, NH white (reference), others (incl. other Hispanic)",
                "education": "3 groups: under high school, high school, above high school (Table 1 'No record' kept missing)",
                "PIR": "< 1, 1 to 3, > 3 (as written)",
                "smoking": "yes = at least 100 cigarettes in life (SMQ020 = 1); fewer = nonsmoker",
                "hypertension": "BPQ020 = 1, mean SBP >= 130 or DBP >= 80, or BPQ050A = 1",
                "hypercholesterolemia": "not defined in the paper: told high cholesterol BPQ080 = 1 (Table 1 has a 'No record' group, i.e. a questionnaire item)",
                "diabetes": "DIQ010 = 1, HbA1c > 6.5 (strict, as written), random serum glucose (LBXSGL) >= 200 mg/dL (11.1 mmol/L), pills (DIQ070; DID070 in 2005-2008) or insulin (DIQ050)",
                "CVD": "cvd_history (MCQ160B-F: CHF, CHD, angina, heart attack, stroke)",
                "BUN, ALT, AST": "continuous as published (LBXSBU mg/dL, LBXSATSI, LBXSASSI U/L)"},
            "notes": "weighted per spec; quartile cut points in the paper not used (continuous target)"}

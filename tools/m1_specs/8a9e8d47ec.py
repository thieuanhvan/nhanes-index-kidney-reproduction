"""M1 builder, 8a9e8d47ec (source withheld). WWI and albuminuria in hypertension.
Spec: WWI per unit (reading_1: WWI is the main exposure), albuminuria = UACR >= 30 mg/g, Model II (Table 2),
NHANES 2005-2006 to 2017-March 2020, no weighting reported.
Sources in fulltext_txt/8a9e8d47ec.txt: 'Definitions and variables' (line 20), Table 1 (variable coding), Table 2 footnote.
"""
import numpy as np
import pandas as pd

REC_ID = "8a9e8d47ec"
CYCLES = ["2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2020"]


def _drink3(lib, df):
    """never / former / current. To 2016: never = < 12 drinks in lifetime (ALQ101 or ALQ100/ALD100 = 2 and ALQ110 = 2);
    ever = >= 12 drinks (any year or lifetime); current = ever and ALQ120Q > 0 (drank in past 12 months), former = ever
    and ALQ120Q = 0. 2017-2020 (no 12-drink items): never = ALQ111 2 (never had a drink), former = ALQ111 1 and ALQ121 0,
    current = ALQ121 1-10."""
    S = lib.stems("ALQ")
    a101 = lib.raw(df, S, ["ALQ101", "ALQ100", "ALD100"])
    a110 = lib.raw(df, S, "ALQ110")
    q120 = lib.raw(df, S, "ALQ120Q").where(lambda v: v <= 365)
    ever = (a101 == 1) | (a110 == 1)
    old = np.select([(a101 == 2) & (a110 == 2), ever & (q120 > 0), ever & (q120 == 0)],
                    ["never", "current", "former"], default="")
    a111 = lib.raw(df, S, "ALQ111")
    a121 = lib.raw(df, S, "ALQ121")
    new = np.select([a111 == 2, a121.between(1, 10), (a111 == 1) & (a121 == 0)], ["never", "current", "former"], default="")
    isnew = df.cycle.isin(["2017-2018", "2017-2020"]).to_numpy()
    v = np.where(isnew, new, old)
    return pd.Series(v, index=df.index).replace("", np.nan).astype(object)


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["WWI"][0](df)
    df["y"] = (df.uacr >= 30).astype(float).where(df.uacr.notna())
    df["drink3"] = _drink3(lib, df)
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["sbp", "dbp", "bmi", "hdl", "ldl", "tc", "crp_mgL", "hba1c", "glucose", "scr"]
    terms += lib.dummies(df, "smoking", "never")
    terms += lib.dummies(df, "drink3", "never")
    # hypertension: self-reported history, antihypertensive medication, SBP >= 140 or DBP >= 90 = hypertension_any_140
    dom = (df.age >= 20) & (df.hypertension_any_140 == 1) & df.x.notna() & df.y.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": None,
            "readings": {
                "which index is the main exposure": "reading_1: WWI, per unit",
                "survey weighting": "none reported (spec weight_variable null): unweighted",
                "age limit": "not stated in the text (Figure 1 only): adults >= 20, the usual NHANES adult cut",
                "pregnancy": "no exclusion stated: not excluded",
                "hypertension (population)": "hypertension_any_140: BPQ020 1 or BPQ050A 1 or mean SBP >= 140 or DBP >= 90",
                "albuminuria": "UACR >= 30 mg/g",
                "race": "RIDRETH1 5 groups (Table 1), ref non-Hispanic White",
                "smoke": "never/former/current from SMQ020, SMQ040",
                "alcohol": "never/former/current ('classified into the same categories'), see _drink3",
                "continuous covariates": "age, SBP, DBP, BMI, HDL-c, LDL-c, TC, CRP, HbA1c, FPG, Cr entered linearly; units "
                                         "of NHANES (mg/dL, mg/L) differ from the paper's mmol/L, umol/L by a constant "
                                         "only, which does not change the WWI coefficient",
                "FPG": "fasting plasma glucose (LBXGLU, fasting subsample) as written",
                "CRP": "CRP (2005-2010) and hs-CRP (2015-2020) in mg/L; not measured 2011-2014"},
            "notes": "the paper does not report missing-data handling; Table 1 shows CRP and FPG summaries for all 15,360 "
                     "(CRP upper quartile 2.63 in every column, suggesting filled values), so the complete-case fit will "
                     "lose 2011-2014 (no CRP) and non-fasting rows (FPG, LDL)."}

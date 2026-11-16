"""M1 builder, e76d3d1271 (Int Urol Nephrol 2024, doi [DOI withheld]). CMI and CKD.
Spec: CMI per 1 unit (reading_1: Table 2 continuous estimate), CKD = eGFR < 60 or albumin >= 30 mg/day (reading_1),
Model 3 (Table 2), NHANES 1999-2000 to 2017-March 2020, NHANES recommended weights.
Sources in fulltext_txt/e76d3d1271.txt: 'Study design and population' (exclusions), 'Calculation of the CMI',
'Primary outcome', 'Definitions of variables of interest', 'Statistical analyses', Table 1 (units), Table 2.
"""
import numpy as np
import pandas as pd

REC_ID = "e76d3d1271"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014",
          "2015-2016", "2017-2020"]
MGDL_RATIO = 88.57 / 38.67   # TG/HDL-C ratio in mg/dL = 2.2904 x the ratio in mmol/L


def _flag(pos, observed):
    return pd.Series(np.where(pos, 1.0, np.where(observed, 0.0, np.nan)), index=pos.index)


def _drink12(lib, df):
    """Yes/no: >= 12 drinks in a year. To 2016: ALQ101 (ALQ100 1999-2000, ALD100 2001-2002) >= 12 drinks in any one
    year; 2017-2020 (no ALQ101): past-year drinks = drinking days (ALQ121 midpoint) x drinks per day (ALQ130) > 12."""
    a101 = lib.raw(df, lib.stems("ALQ"), ["ALQ101", "ALQ100", "ALD100"])
    out = pd.Series(np.where(a101 == 1, 1.0, np.where(a101 == 2, 0.0, np.nan)), index=df.index)
    days = {0: 0, 1: 365, 2: 286, 3: 182, 4: 104, 5: 52, 6: 30, 7: 12, 8: 9, 9: 4.5, 10: 1.5}
    a111 = lib.raw(df, lib.stems("ALQ"), "ALQ111")
    a121 = lib.raw(df, lib.stems("ALQ"), "ALQ121").map(days)
    a130 = lib.raw(df, lib.stems("ALQ"), "ALQ130").where(lambda v: v < 77)
    new = df.cycle.isin(["2017-2018", "2017-2020"])
    annual = a121 * a130.where(a121 > 0, 0.0)
    v = np.where(a111 == 2, 0.0, np.where(annual.notna(), (annual > 12).astype(float), np.nan))
    out[new] = v[new.to_numpy()]
    return out


def build(lib):
    df = lib.table(CYCLES)
    # CMI with TG/HDL-C in mg/dL: Table 1 gives TG and HDL-C in mg/dl and TG/HDL-C mean 2.89, CMI mean 1.67 = 2.89 x WHtR 0.58
    df["x"] = lib.INDEX["CMI"][0](df) * MGDL_RATIO
    ok = df.uacr.notna() & df.egfr_ckdepi2009.notna()
    df["y"] = ((df.uacr >= 30) | (df.egfr_ckdepi2009 < 60)).astype(float).where(ok)
    S = lib.stems
    df["smoke_yes"] = lib.raw(df, S("SMQ"), "SMQ020").map({1: 1.0, 2: 0.0})
    df["drink_yes"] = _drink12(lib, df)
    # hypertension: previous diagnosis, or SBP >= 130 and/or DBP >= 80
    hpos = (df.BPQ020 == 1) | (df.sbp >= 130) | (df.dbp >= 80)
    df["htn"] = _flag(hpos, df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna()))
    # DM: self-reported diabetes history or fasting glucose >= 7.0 mmol/L (126 mg/dL)
    dpos = (df.DIQ010 == 1) | (df.glucose >= 126)
    df["dm"] = _flag(dpos, df.DIQ010.isin([1, 2, 3]) | df.glucose.notna())
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["smoke_yes", "drink_yes", "bmi", "tc", "ldl", "htn", "dm"]
    df["w"] = lib.pooled_weight(df, "fasting")
    cancer = lib.raw(df, S("MCQ"), "MCQ220")
    dom = (df.age > 20) & (df.pregnant != 1) & (cancer != 1) & df.x.notna() & df.y.notna() & (df.w > 0)
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "exposure scale of the abstract estimate": "reading_1: CMI per 1 unit, continuous",
                "albuminuria component of CKD definition": "reading_1 (>= 30 mg albumin per day): NHANES has only spot "
                    "urine, so the only public-file measure is UACR >= 30 mg/g (KDIGO equivalent of 30 mg/day); this "
                    "coincides with reading_2",
                "CKD 'consistently low eGFR'": "single visit eGFR < 60",
                "eGFR equation": "not stated: CKD-EPI 2009 with race",
                "CMI unit": "TG/HDL-C ratio in mg/dL (Table 1: TG, HDL-C in mg/dl, TG/HDL-C 2.89, CMI 1.67); "
                            "lib.INDEX['CMI'] (mmol/L ratio) multiplied by 88.57/38.67",
                "age": "'above 20 years of age' read as age > 20",
                "cancer": "MCQ220 = 1 (ever told cancer or malignancy) excluded",
                "survey weight": "fasting subsample weight (CMI needs fasting TG; LDL-C covariate); weight 0 excluded",
                "race": "RIDRETH1 5 groups (Table 1), ref non-Hispanic White",
                "smoke": "yes = SMQ020 1 (>= 100 cigarettes; Table 1 'Smoke' 45.7 %)",
                "drink": "yes = >= 12 drinks in a year (ALQ101 family; 2017-2020 past-year drinks > 12); definition not "
                         "given in the paper",
                "hypertension": "BPQ020 1 or SBP >= 130 or DBP >= 80 (medication not in the definition)",
                "DM": "DIQ010 1 or FPG >= 126 mg/dL",
                "BMI, TC, LDL-C, age": "continuous"},
            "notes": "CMI unit chosen from Table 1 before the first execution; it differs from the recon key CMI (mmol/L), "
                     "flag for the main session."}

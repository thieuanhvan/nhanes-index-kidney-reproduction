"""M1 builder, 8991f10763 (Medicine 2025, doi [DOI withheld]). AISI and CKD.
Spec: AISI/100 per unit, CKD = eGFR < 60 or ACR >= 30 mg/g, fully adjusted model 3 (Table 2), NHANES 2007-2018, weighted.
Sources in fulltext_txt/8991f10763.txt: sections 2.2 (exclusions), 2.3 (AISI), 2.4 (CKD), 2.5 (covariates), 2.6
(weights, models), Table 2 footnote (model 3 set, spec reading_1: hyperuricemia, no blood uric acid).
"""
import numpy as np
import pandas as pd

REC_ID = "8991f10763"
CYCLES = ["2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def _cat(df, name, conds, labels, base):
    df[name] = pd.Series(np.select(conds, labels, default=""), index=df.index).astype(object)
    df.loc[df[base].isna(), name] = np.nan
    return name


def _flag(pos, observed):
    return pd.Series(np.where(pos, 1.0, np.where(observed, 0.0, np.nan)), index=pos.index)


def _drink12(lib, df):
    """Yes/no: more than 12 drinks a year. 2007-2016 ALQ101 (>= 12 drinks in any one year); 2017-2018 (no ALQ101):
    past-year drinks = drinking days (ALQ121 midpoint) x drinks per day (ALQ130) > 12, ALQ111 2 = no."""
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
    df["x"] = lib.INDEX["AISI"][0](df) / 100
    ok = df.uacr.notna() & df.egfr_ckdepi2009.notna()
    df["y"] = ((df.uacr >= 30) | (df.egfr_ckdepi2009 < 60)).astype(float).where(ok)
    S = lib.stems
    _cat(df, "age3", [df.age < 40, df.age < 60, df.age >= 60], ["20_39", "40_59", "ge60"], "age")
    _cat(df, "pir3", [df.pir < 1.3, df.pir <= 3.5, df.pir > 3.5], ["lt1.3", "1.3_3.5", "gt3.5"], "pir")
    _cat(df, "bmi3", [df.bmi < 25, df.bmi <= 30, df.bmi > 30], ["lt25", "25_30", "gt30"], "bmi")
    df["marital2"] = lib.raw(df, S("DEMO"), "DMDMARTL").map({1: "coupled", 6: "coupled", 2: "single_sep",
                                                            3: "single_sep", 4: "single_sep", 5: "single_sep"})
    df["smoke_yes"] = lib.raw(df, S("SMQ"), "SMQ020").map({1: 1.0, 2: 0.0})
    df["drink_yes"] = _drink12(lib, df)
    # hyperuricemia: serum urate > 7.0 mg/dL men, > 5.7 mg/dL women
    df["hyperuric"] = np.where(df.female == 1, df.uric_acid > 5.7, df.uric_acid > 7.0).astype(float)
    df.loc[df.uric_acid.isna() | df.female.isna(), "hyperuric"] = np.nan
    # hypertension: mean SBP >= 140 or mean DBP >= 90, or self-report along with antihypertensive drugs
    b050a = lib.raw(df, S("BPQ"), "BPQ050A")
    hpos = (df.sbp >= 140) | (df.dbp >= 90) | ((df.BPQ020 == 1) & (b050a == 1))
    df["htn"] = _flag(hpos, (df.sbp.notna() & df.dbp.notna()) | df.BPQ020.isin([1, 2]))
    # DM: HbA1c > 6.5, FPG >= 7.0 mmol/L (126 mg/dL), diagnosis, DM medication or insulin
    pills = lib.raw(df, S("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, S("DIQ"), "DIQ050")
    dpos = (df.hba1c > 6.5) | (df.glucose >= 126) | (df.DIQ010 == 1) | (pills == 1) | (ins == 1)
    df["dm"] = _flag(dpos, df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna())
    terms = ["x"]
    terms += lib.dummies(df, "age3", "20_39")
    terms += ["female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "marital2", "coupled")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "pir3", "lt1.3")
    terms += lib.dummies(df, "bmi3", "lt25")
    terms += ["alt", "ast", "tc", "tg_serum", "drink_yes", "smoke_yes", "htn", "hyperuric", "dm"]
    df["w"] = lib.pooled_weight(df, "mec")
    dom = (df.age >= 20) & df.y.notna() & df.x.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "Model 3 covariate set": "reading_1: Table 2 footnote, hyperuricemia, no blood uric acid",
                "eGFR equation": "CKD-EPI 2009 with race (version not stated)",
                "ACR threshold": ">= 30 mg/g as written",
                "survey weight": "MEC weight (AISI from CBC; variable not named)",
                "age": "categories 20-39, 40-59, >= 60 as section 2.5 codes it",
                "PIR": "< 1.3, 1.3-3.5 (inclusive), > 3.5",
                "BMI": "< 25, 25-30.0 (inclusive), > 30",
                "marital": "coupled = married or living with partner (DMDMARTL 1, 6); single/separated = 2-5",
                "education": "educ3: below HS (1-2), HS/GED (3), above HS (4-5)",
                "race": "RIDRETH1 5 groups, ref non-Hispanic White",
                "smoking": "yes = SMQ020 1 (> 100 cigarettes lifetime)",
                "drinking": "yes = ALQ101 1 2007-2016; 2017-2018 past-year drinks (ALQ121 day midpoint x ALQ130) > 12",
                "TG": "serum TG from the biochemistry profile (LBXSTR, all MEC examinees): the paper does not say fasting "
                      "and analyses 24,386 of 34,770 adults, more than the fasting subsample holds",
                "TC, AST, ALT": "continuous as reported (mg/dL, IU/L)",
                "hypertension": "SBP >= 140 or DBP >= 90 or (BPQ020 1 and BPQ050A 1): 'self-report along with the use "
                                "of antihypertensive drugs' read literally as both",
                "hyperuricemia": "> 7.0 mg/dL men, > 5.7 mg/dL women ('exceeding')",
                "DM": "HbA1c > 6.5 ('exceeding'), FPG >= 126 mg/dL, DIQ010 1, pills (DIQ070, DID070 2007-2008) or "
                      "insulin (DIQ050)"},
            "notes": "domain = adults >= 20 with CKD and AISI data; covariate completeness applied by the fit."}

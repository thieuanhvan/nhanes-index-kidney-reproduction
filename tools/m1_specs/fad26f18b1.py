"""M1 builder, fad26f18b1 (PLOS One 2026, doi [DOI withheld]). HALP score and CKD.
Spec: HALP per 1 point, CKD = ACR >= 30 or eGFR (CKD-EPI 2009) < 60, Model 3 (Table 2), NHANES 2009-2018, weighted.
Sources in fulltext_txt/fad26f18b1.txt: exclusions lines 145-151, HALP 152-156, CKD 157-162, covariates 167-191,
statistics 252-268 (multiple imputation of missing covariates; complete case here, applied by the fit).
"""
import numpy as np
import pandas as pd

REC_ID = "fad26f18b1"
CYCLES = ["2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def _cat(df, name, conds, labels, base):
    df[name] = pd.Series(np.select(conds, labels, default=""), index=df.index).astype(object)
    df.loc[df[base].isna(), name] = np.nan
    return name


def _drink12(lib, df):
    """Yes/no: at least 12 alcoholic drinks in a year. 2009-2016 ALQ101 (had >= 12 drinks in any one year);
    2017-2018 (no ALQ101): past-year drinks = drinking days (ALQ121 midpoint) x drinks per day (ALQ130) > 12."""
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


def _flag(pos, observed):
    return pd.Series(np.where(pos, 1.0, np.where(observed, 0.0, np.nan)), index=pos.index)


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["HALP"][0](df)
    ok = df.uacr.notna() & df.egfr_ckdepi2009.notna()
    df["y"] = ((df.uacr >= 30) | (df.egfr_ckdepi2009 < 60)).astype(float).where(ok)
    S = lib.stems
    # demographic covariates as the Covariates section codes them (spec reading_1)
    _cat(df, "age2", [df.age < 60, df.age >= 60], ["lt60", "ge60"], "age")
    _cat(df, "pir3", [df.pir < 1, df.pir <= 4, df.pir > 4], ["low", "mid", "high"], "pir")
    _cat(df, "bmi3", [df.bmi < 25, df.bmi < 30, df.bmi >= 30], ["lt25", "25_30", "ge30"], "bmi")
    df["educ2"] = df.educ_raw.map({1: "hs_or_less", 2: "hs_or_less", 3: "hs_or_less", 4: "gt_hs", 5: "gt_hs"})
    df["marital"] = lib.raw(df, S("DEMO"), "DMDMARTL").map({1: "married", 2: "widowed", 3: "divorced", 4: "separated",
                                                           5: "never", 6: "partner"})
    smq = lib.raw(df, S("SMQ"), "SMQ020")
    df["smoke_yes"] = smq.map({1: 1.0, 2: 0.0})
    df["drink_yes"] = _drink12(lib, df)
    # comorbidities as defined in the paper
    # DM: HbA1c >= 6.5, FPG >= 7.0 mmol/L (126 mg/dL), diagnosis, DM medication or insulin (pills DIQ070, DID070 in 2005-2008)
    pills = lib.raw(df, S("DIQ"), ["DIQ070", "DID070"])
    ins = lib.raw(df, S("DIQ"), "DIQ050")
    dpos = (df.hba1c >= 6.5) | (df.glucose >= 126) | (df.DIQ010 == 1) | (pills == 1) | (ins == 1)
    df["dm"] = _flag(dpos, df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna())
    # hypertension: mean SBP >= 140, mean DBP >= 90, self-reported, or current antihypertensive use = hypertension_any_140
    # CVD: congestive heart failure, coronary heart disease, angina or heart attack (MCQ160B-E; stroke not listed)
    mc = df[["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E"]]
    df["cvd4"] = _flag((mc == 1).any(axis=1), mc.isin([1, 2]).all(axis=1))
    # dyslipidemia: TC >= 5.2, TG >= 1.7, LDL >= 3.4, HDL <= 1.0 mmol/L, self-reported high cholesterol, lipid-lowering drug
    b080 = lib.raw(df, S("BPQ"), "BPQ080")
    b100d = lib.raw(df, S("BPQ"), "BPQ100D")
    lpos = ((df.tc >= 5.2 * 38.67) | (df.tg >= 1.7 * 88.57) | (df.ldl >= 3.4 * 38.67) | (df.hdl <= 1.0 * 38.67)
            | (b080 == 1) | (b100d == 1))
    df["dyslip"] = _flag(lpos, (df.tc.notna() & df.hdl.notna()) | b080.isin([1, 2]))
    terms = ["x"]
    terms += lib.dummies(df, "age2", "lt60")
    terms += ["female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "marital", "married")
    terms += lib.dummies(df, "educ2", "hs_or_less")
    terms += lib.dummies(df, "pir3", "low")
    terms += lib.dummies(df, "bmi3", "lt25")
    terms += ["smoke_yes", "drink_yes", "ast", "alt", "hypertension_any_140", "dm", "cvd4", "dyslip"]
    df["w"] = lib.pooled_weight(df, "mec")
    dial = lib.raw(df, S("KIQ_U"), "KIQ025")
    dom = ((df.age >= 20) & (df.pregnant != 1) & df.scr.notna() & df.uacr.notna() & df.x.notna() & (dial != 1))
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "Demographic covariate coding": "reading_1: age <60/>=60, PIR <1, 1-4, >4, BMI <25, 25-<30, >=30",
                "CKD definition by single visit": "reading_1: ACR >= 30 or eGFR < 60 at one visit",
                "eGFR equation": "CKD-EPI 2009 with race (stated 2009)",
                "survey weight": "MEC weight (HALP from CBC and biochemistry; variable not named)",
                "race": "RIDRETH1 5 groups, ref non-Hispanic White",
                "education": "high school or less (DMDEDUC2 1-3) vs more (4-5)",
                "marital": "DMDMARTL 6 categories, ref married",
                "smoking": "yes = SMQ020 1 (>= 100 cigarettes lifetime)",
                "drinking": "yes = ALQ101 1 (>= 12 drinks in any one year) 2009-2016; 2017-2018 past-year drinks "
                            "(ALQ121 day midpoint x ALQ130) > 12, ALQ111 2 = no",
                "DM": "HbA1c >= 6.5 or FPG >= 126 mg/dL or DIQ010 1 or pills (DIQ070) or insulin (DIQ050)",
                "hypertension": "hypertension_any_140 (BPQ020, BPQ050A, SBP >= 140, DBP >= 90)",
                "CVD": "MCQ160B-E any yes (CHF, CHD, angina, heart attack); stroke not in the paper's list",
                "dyslipidemia": "TC >= 201.0, TG >= 150.5, LDL >= 131.4, HDL <= 38.6 mg/dL (mmol cut points x factors), "
                                "BPQ080 1 or BPQ100D 1; negative when TC and HDL observed or BPQ080 answered",
                "dialysis": "KIQ025 = 1 (dialysis in past 12 months) excluded",
                "AST, ALT": "continuous U/L"},
            "notes": "paper used multiple imputation for missing covariates; M1 fit is complete case. "
                     "TG and LDL are fasting-subsample variables, used only inside the dyslipidemia flag."}

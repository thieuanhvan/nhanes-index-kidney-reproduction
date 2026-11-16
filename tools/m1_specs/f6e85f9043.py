"""M1 builder, f6e85f9043 (CTI and albuminuria, Medicine 2025, doi [DOI withheld]).
Spec: CTI per 1 unit, albuminuria UACR >= 30 mg/g, Model 3 (Table 2 footnote, 18 covariates), NHANES 1999-2010, no weighting.
Sources in the paper (fulltext_txt/f6e85f9043.txt): exclusions section 2.1 (line 12), CTI 2.2 (line 15), outcome 2.3 (line 18),
covariates 2.4 (line 21), models 2.5 (line 25), Table 1 categories (lines 105-163), Table 2 footnote (line 178).
"""
import numpy as np
import pandas as pd

REC_ID = "f6e85f9043"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010"]


def _miss(v, codes=(7, 9, 77, 99)):
    return v.where(~v.isin(codes))


def _cat(df, name, conds, labels, known):
    df[name] = np.select(conds, labels, default="").astype(object)
    df.loc[~known, name] = np.nan
    df.loc[df[name] == "", name] = np.nan


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    # exposure: CTI = 0.412 ln(CRP mg/dL) + ln(TG mg/dL * FPG mg/dL / 2)
    df["x"] = lib.INDEX["CTI_mgdl"][0](df).replace([np.inf, -np.inf], np.nan)
    df["y"] = (df.uacr >= 30).astype(float).where(df.uacr.notna())

    # age (<60, >=60) as in section 2.4
    df["age_ge60"] = (df.age >= 60).astype(float).where(df.age.notna())
    # marital: married or living with partner vs living alone (Table 1)
    mar = _miss(lib.raw(df, s("DEMO"), "DMDMARTL"))
    df["married"] = mar.isin([1, 6]).astype(float).where(mar.isin([1, 2, 3, 4, 5, 6]))
    # PIR categories <1.3, 1.3-3.5, >=3.5 (Table 1)
    _cat(df, "pir_cat", [df.pir < 1.3, df.pir < 3.5, df.pir >= 3.5], ["lt1.3", "1.3-3.5", "ge3.5"], df.pir.notna())
    # BMI categories <25, 25-30, >=30 (Table 1)
    _cat(df, "bmi_cat", [df.bmi < 25, df.bmi < 30, df.bmi >= 30], ["lt25", "25-30", "ge30"], df.bmi.notna())
    # low eGFR: eGFR < 60, CKD-EPI (version not stated: 2009 with race)
    df["low_egfr"] = (df.egfr_ckdepi2009 < 60).astype(float).where(df.egfr_ckdepi2009.notna())
    # diabetes: physician diagnosis, FPG >= 7.0 mmol/L, HbA1c >= 6.5 %, OGTT 2 h >= 11.1 mmol/L (OGTT files 2005-2010)
    ogtt = lib.raw(df, {"2005-2006": "OGTT_D", "2007-2008": "OGTT_E", "2009-2010": "OGTT_F"}, "LBXGLT")
    pos = (df.DIQ010 == 1) | (df.glucose >= 126) | (df.hba1c >= 6.5) | (ogtt >= 200)
    df["dm"] = pd.Series(np.where(pos, 1.0, np.where(df.DIQ010.isin([1, 2, 3]) | df.glucose.notna() | df.hba1c.notna(), 0.0, np.nan)), index=df.index)
    # alcohol: nondrinker (<12 drinks lifetime and in any one year), former (>=12, none past 12 months), current
    a12y = _miss(lib.raw(df, s("ALQ"), ["ALQ101", "ALD100", "ALQ100"]))
    a12l = _miss(lib.raw(df, s("ALQ"), "ALQ110"))
    freq = lib.raw(df, s("ALQ"), "ALQ120Q")
    freq = freq.where(~freq.isin([777, 999]))
    ever = (a12y == 1) | (a12l == 1)
    never = (a12y == 2) & (a12l == 2)
    _cat(df, "alcohol", [never, ever & (freq == 0), ever & (freq > 0)], ["never", "former", "current"],
         never | (ever & freq.notna()))
    # physical activity: vigorous (yes/no), moderate (yes/no)
    # 1999-2006 PAD200 / PAD320 (leisure, past 30 days; 3 'unable' = no); 2007-2010 GPAQ work or recreation
    pad200 = _miss(lib.raw(df, s("PAQ"), "PAD200"))
    pad320 = _miss(lib.raw(df, s("PAQ"), "PAD320"))
    p605, p650 = _miss(lib.raw(df, s("PAQ"), "PAQ605")), _miss(lib.raw(df, s("PAQ"), "PAQ650"))
    p620, p665 = _miss(lib.raw(df, s("PAQ"), "PAQ620")), _miss(lib.raw(df, s("PAQ"), "PAQ665"))
    old = df.cycle.isin(CYCLES[:4])
    vig_new = ((p605 == 1) | (p650 == 1)).astype(float).where(p605.notna() & p650.notna() | (p605 == 1) | (p650 == 1))
    mod_new = ((p620 == 1) | (p665 == 1)).astype(float).where(p620.notna() & p665.notna() | (p620 == 1) | (p665 == 1))
    df["pa_vig"] = np.where(old, (pad200 == 1).astype(float).where(pad200.isin([1, 2, 3])), vig_new)
    df["pa_mod"] = np.where(old, (pad320 == 1).astype(float).where(pad320.isin([1, 2, 3])), mod_new)

    terms = ["x", "age_ge60", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "bmi_cat", "lt25")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["married"]
    terms += lib.dummies(df, "pir_cat", "lt1.3")
    terms += ["low_egfr", "uric_acid", "albumin_serum", "tg", "ldl", "dm"]
    terms += lib.dummies(df, "alcohol", "never")
    terms += ["pa_vig", "pa_mod"]
    terms += lib.dummies(df, "smoking", "never")
    terms += ["alt", "ast"]

    dom = (df.age >= 20) & (df.pregnant != 1) & df.x.notna() & df.uacr.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": None,
            "readings": {
                "covariate list": "Table 2 footnote Model 3, 18 covariates, no hypertension (reading_1)",
                "albuminuria threshold": "UACR >= 30 mg/g (as written)",
                "age limit": "none in Methods; abstract 'US adults' read as age >= 20 (NHANES adult convention, recon CFG)",
                "survey weighting": "none stated (weight_variable null): unweighted",
                "eGFR equation (low eGFR)": "CKD-EPI 2009 with race (version not stated)",
                "age coding": "binary <60 / >=60 as stated in section 2.4",
                "BMI and PIR coding": "categories of Table 1 (BMI <25, 25-30, >=30; PIR <1.3, 1.3-3.5, >=3.5)",
                "education": "3 levels (under high school / high school / college or above), DMDEDUC2 1-2, 3, 4-5",
                "marital status": "married or living with partner (DMDMARTL 1, 6) vs living alone (2-5)",
                "race": "RIDRETH1 5 groups (Table 1)",
                "alcohol": "never (ALQ101/ALD100/ALQ100 = 2 and ALQ110 = 2), former (>=12 drinks, ALQ120Q = 0), current (ALQ120Q > 0)",
                "physical activity": "two binary terms vigorous yes/no and moderate yes/no; 1999-2006 PAD200/PAD320 (unable = no), 2007-2010 GPAQ work or recreation",
                "diabetes": "DIQ010 = 1, FPG >= 126 mg/dL, HbA1c >= 6.5 %, OGTT 2 h >= 200 mg/dL (2005-2010 only); insulin/pills not in the paper's definition",
                "TG, LDL-C": "fasting LBXTR / LBDLDL as continuous mg/dL",
                "missing data": "paper used multiple imputation; M1 fit is complete case (automatic)",
            },
            "notes": "builder written before first execution; CTI from recon INDEX CTI_mgdl (CRP mg/L / 10)"}

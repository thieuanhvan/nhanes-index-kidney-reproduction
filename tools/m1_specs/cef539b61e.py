"""M1 builder, cef539b61e (SII and CKD, Front Immunol 2024, doi [DOI withheld]).
Spec: SII/1000 per 1 unit, CKD (eGFR < 60 or UACR >= 30), Model 3 (24 covariates), NHANES 2007-2018, weighted.
Sources in the paper (fulltext_txt/cef539b61e.txt): population line 14, exposure line 19, outcome line 22,
covariates line 25, models line 28, Table 1 categories lines 78-131, Table 2 reference levels lines 133-170.
"""
import numpy as np
import pandas as pd

REC_ID = "cef539b61e"
CYCLES = ["2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]
OGTT = {"2007-2008": "OGTT_E", "2009-2010": "OGTT_F", "2011-2012": "OGTT_G", "2013-2014": "OGTT_H", "2015-2016": "OGTT_I"}


def _miss(v, codes=(7, 9, 77, 99)):
    return v.where(~v.isin(codes))


def _alcohol5(df, lib):
    """never / former / mild / moderate / heavy (usual NHANES grouping, e.g. Rattan et al.):
    never: < 12 drinks in lifetime (2007-2016 ALQ101 = 2 and ALQ110 = 2; 2017-2018 ALQ111 = 2, never had a drink)
    former: >= 12 drinks (ever drank in 2017-2018) but none in the past 12 months
    heavy: >= 3 drinks/day women, >= 4 men, or binge on >= 5 days per month
    moderate: >= 2 drinks/day women, >= 3 men, or binge on >= 2 days per month
    mild: other current drinkers."""
    s = lib.stems("ALQ")
    a101 = _miss(lib.raw(df, s, "ALQ101"))
    a110 = _miss(lib.raw(df, s, "ALQ110"))
    q120 = lib.raw(df, s, "ALQ120Q"); q120 = q120.where(~q120.isin([777, 999]))
    u120 = _miss(lib.raw(df, s, "ALQ120U"))
    dpd = lib.raw(df, s, "ALQ130"); dpd = dpd.where(~dpd.isin([777, 999]))
    # binge days: 2007-2010 ALQ140Q/U (5+ drinks), 2011-2016 ALQ141Q/U (4/5+ drinks); unit 1 week, 2 month, 3 year
    bq = lib.raw(df, s, ["ALQ141Q", "ALQ140Q"]); bq = bq.where(~bq.isin([777, 999]))
    bu = _miss(lib.raw(df, s, ["ALQ141U", "ALQ140U"]))
    per_month = {1: 52 / 12, 2: 1.0, 3: 1 / 12}
    binge = bq * bu.map(per_month)
    binge = binge.where(~(bq == 0), 0.0)
    # 2017-2018 items
    a111 = _miss(lib.raw(df, s, "ALQ111"))
    a121 = _miss(lib.raw(df, s, "ALQ121"))
    a142 = _miss(lib.raw(df, s, "ALQ142"))
    code_month = {0: 0.0, 1: 30.0, 2: 25.0, 3: 15.0, 4: 8.6, 5: 4.3, 6: 2.5, 7: 1.0, 8: 0.75, 9: 0.375, 10: 0.125}
    new = df.cycle == "2017-2018"
    never = np.where(new, a111 == 2, (a101 == 2) & (a110 == 2))
    ever = np.where(new, a111 == 1, (a101 == 1) | (a110 == 1))
    cur_freq = np.where(new, a121, q120)
    binge_m = np.where(new, a142.map(code_month), binge)
    cur_freq = pd.Series(cur_freq, index=df.index)
    binge_m = pd.Series(binge_m, index=df.index)
    fem = df.female == 1
    current = ever & (cur_freq > 0)
    heavy = current & (((dpd >= 3) & fem) | ((dpd >= 4) & ~fem) | (binge_m >= 5))
    moder = current & ~heavy & (((dpd >= 2) & fem) | ((dpd >= 3) & ~fem) | (binge_m >= 2))
    mild = current & ~heavy & ~moder & dpd.notna() & binge_m.notna()
    former = ever & (cur_freq == 0)
    out = pd.Series(np.select([never, former, heavy, moder, mild], ["never", "former", "heavy", "moderate", "mild"],
                              default=""), index=df.index).astype(object)
    return out.where(out != "", np.nan)


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["SII"][0](df).replace([np.inf, -np.inf], np.nan) / 1000
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())

    df["alcohol5"] = _alcohol5(df, lib)
    # hyperlipidemia (not defined in the paper; NCEP ATP III convention): TG >= 150, TC >= 200, LDL >= 130,
    # HDL < 40 men / < 50 women, or taking prescribed cholesterol medicine (BPQ100D = 1)
    bpq100d = _miss(lib.raw(df, lib.stems("BPQ"), "BPQ100D"))
    hdl_low = np.where(df.female == 1, df.hdl < 50, df.hdl < 40)
    pos = (df.tg >= 150) | (df.tc >= 200) | (df.ldl >= 130) | hdl_low | (bpq100d == 1)
    df["hyperlip"] = pd.Series(np.where(pos, 1.0, np.where(df.tc.notna() & df.hdl.notna(), 0.0, np.nan)), index=df.index)
    # diabetes 3 levels (Table 1): DM = diagnosis, insulin, pills, HbA1c >= 6.5, FPG >= 126, OGTT >= 200;
    # pre-DM = told prediabetes (DIQ160), FPG 100-125, HbA1c 5.7-6.4, OGTT 140-199; else no
    ogtt = lib.raw(df, OGTT, "LBXGLT")
    diq160 = _miss(lib.raw(df, lib.stems("DIQ"), "DIQ160"))
    dm = (df.diabetes_any == 1) | (ogtt >= 200)
    pre = (diq160 == 1) | df.glucose.between(100, 125.999) | df.hba1c.between(5.7, 6.49) | ogtt.between(140, 199.999)
    known = df.diabetes_any.notna() | ogtt.notna()
    df["dm3"] = pd.Series(np.select([dm, pre], ["dm", "predm"], default="no"), index=df.index).astype(object).where(known, np.nan)

    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "alcohol5", "never")
    terms += lib.dummies(df, "smoking", "never")
    terms += ["pir", "bmi", "waist", "dbp", "sbp", "glucose", "alt", "ast", "scr", "uric_acid", "tg", "tc", "hdl", "ldl",
              "egfr_ckdepi2009", "hyperlip", "hypertension_any_140"]
    terms += lib.dummies(df, "dm3", "no")

    df["w"] = lib.pooled_weight(df, "fasting")
    dom = (df.age >= 18) & (df.pregnant != 1) & df.x.notna() & df.uacr.notna() & e.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "weight variable": "reading_1: subsample weight; the 17,651 dropped for missing weight and the use of Glu, TG, LDL-C point to the fasting subsample: pooled WTSAF2YR",
                "eGFR equation": "reading_1 CKD-EPI version not stated: 2009 with race (brief rule 5)",
                "CKD threshold": "UACR >= 30 or eGFR < 60 (as written)",
                "exposure scale": "SII/1000 per 1 unit (SII = neutrophils x platelets / lymphocytes, 10^3 cells/uL)",
                "Glu": "fasting plasma glucose LBXGLU (fasting subsample, consistent with the weight)",
                "TG, LDL-C": "fasting LBXTR / LBDLDL; TC, HDL-C from TCHOL / HDL files",
                "race": "RIDRETH1 5 groups (Table 1)",
                "education": "3 levels (Table 1), DMDEDUC2 1-2 / 3 / 4-5",
                "alcohol status": "5 levels of Table 1 (never/former/mild/moderate/heavy), not defined in the paper; usual NHANES grouping (Rattan) coded in _alcohol5",
                "smoking": "never / former / now (nhanes_core smoking)",
                "hypertension": "not defined in the paper; convention: diagnosis, current medication, or mean BP >= 140/90 (hypertension_any_140)",
                "hyperlipidemia": "not defined in the paper; NCEP ATP III convention with BPQ100D for medication",
                "diabetes mellitus": "3 levels no / pre-DM / DM as in Table 1; definitions are the usual NHANES ones (paper gives none)",
                "continuous covariates": "PIR, BMI, WC, DBP, SBP, Glu, ALT, AST, Cr, UA, TG, TC, HDL-C, LDL-C, eGFR entered linearly (units of the public files)",
            },
            "notes": "builder written before first execution"}

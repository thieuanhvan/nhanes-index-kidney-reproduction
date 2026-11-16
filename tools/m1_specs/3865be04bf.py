"""M1 builder, 3865be04bf (FIB4 and CKD, PLoS One 2025, doi [DOI withheld]).
Spec: FIB4 (AST, ALT in ukat/L) per 1 unit, CKD = eGFR < 60 (CKD-EPI 2009) or UACR >= 30, Model 3 (Table 2 footnote, no
marital status, reading_1), NHANES 2005-2006 to 2015-2016 plus 2017-March 2020 (reading_1), weighted.
Sources in fulltext_txt/3865be04bf.txt: population line 12, CKD and FIB4 line 17, eGFR equation lines 21-22, covariates
line 25, weighting and no age adjustment line 28, Table 1 categories lines 141-195, Table 2 footnote line 206.
"""
import numpy as np
import pandas as pd

REC_ID = "3865be04bf"
CYCLES = ["2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2020"]
ALQ121_DAYS = {0: 0, 1: 365, 2: 286, 3: 182, 4: 104, 5: 52, 6: 30, 7: 12, 8: 9, 9: 4.5, 10: 1.5}


def _miss(v, codes):
    return v.where(~v.isin(codes))


def _drink_days(lib, df):
    s = lib.stems("ALQ")
    q = _miss(lib.raw(df, s, "ALQ120Q"), [777, 999])
    u = _miss(lib.raw(df, s, "ALQ120U"), [7, 9])
    days = (q * u.map({1: 52, 2: 12, 3: 1})).where(q != 0, 0.0)
    q121 = _miss(lib.raw(df, s, "ALQ121"), [77, 99])
    days = days.fillna(q121.map(ALQ121_DAYS))
    skip = (lib.raw(df, s, "ALQ110") == 2) | (lib.raw(df, s, "ALQ111") == 2)
    return days.where(~(skip & days.isna()), 0.0)


def _met_week(lib, df):
    """MET-min/week. GPAQ cycles: 8 x vigorous + 4 x moderate (work, leisure) + 4 x walk/bike, missing parts 0.
    2005-2006: leisure activities file PAQIAF, sum of PADMETS x PADTIMES (30 days) x PADDURAT x 7/30."""
    s = lib.stems("PAQ")
    tot = pd.Series(0.0, index=df.index)
    for gate, d, m, met in [("PAQ605", "PAQ610", "PAD615", 8), ("PAQ620", "PAQ625", "PAD630", 4), ("PAQ635", "PAQ640", "PAD645", 4),
                            ("PAQ650", "PAQ655", "PAD660", 8), ("PAQ665", "PAQ670", "PAD675", 4)]:
        g = lib.raw(df, s, gate)
        dd = _miss(lib.raw(df, s, d), [77, 99])
        mm = _miss(lib.raw(df, s, m), [7777, 9999])
        tot += (met * dd * mm).where(g == 1).fillna(0.0)
    m = (df.cycle == "2005-2006").to_numpy()
    x = lib._read("2005-2006", "PAQIAF_D")
    if m.any() and x is not None:
        t = _miss(x["PADTIMES"], [77777, 99999])
        du = _miss(x["PADDURAT"], [77777, 99999])
        wk = (x["PADMETS"] * t * du * 7 / 30).groupby(x["SEQN"]).sum()
        tot[m] = np.nan_to_num(wk.reindex(df.loc[m, "SEQN"]).to_numpy(float))
    return tot.where(df.has_paq.astype(bool))


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["FIB4_ukat"][0](df)
    df["y"] = ((df.egfr_ckdepi2009 < 60) | (df.uacr >= 30)).astype(float).where(df.uacr.notna() & df.egfr_ckdepi2009.notna())
    s = lib.stems
    # liver exclusions definable in the downloaded public files
    mcq170l = lib.raw(df, s("MCQ"), "MCQ170L") == 1
    cirr = lib.raw(df, s("MCQ"), "MCQ510C")
    hep = lib.raw(df, s("MCQ"), "MCQ510D")
    liver = mcq170l | (cirr.notna() & (cirr < 77)) | (hep.notna() & (hep < 77))
    df["poor"] = (df.pir < 1).astype(float).where(df.pir.notna())
    df["ever_smoker"] = (df.SMQ020 == 1).astype(float).where(df.SMQ020.isin([1, 2]))
    dd = _drink_days(lib, df)
    df["drink_cur"] = (dd > 0).astype(float).where(dd.notna())
    met = _met_week(lib, df)
    df["low_pa"] = (met < 600).astype(float).where(met.notna())
    pills = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    ogtt = lib.raw(df, s("OGTT"), "LBXGLT")
    dpos = ((df.DIQ010 == 1) | (df.DIQ050 == 1) | (pills == 1) | (df.glucose >= 126) | (df.glucose_serum >= 200)
            | (ogtt >= 200))
    df["diab"] = dpos.astype(float).where(dpos | df.DIQ010.isin([1, 2, 3]))
    terms = ["x", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["poor", "ever_smoker", "drink_cur", "low_pa"]
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["hypertension_any_140", "diab"]
    df["w"] = lib.pooled_weight(df, "mec")
    dom = (df.age >= 18) & ~liver & df.x.notna() & df.y.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "marital status": "not in Model 3 (Table 2 footnote, reading_1)",
                "cycles": "2005-2006 to 2015-2016 plus 2017-March 2020 P file (reading_1)",
                "eGFR": "CKD-EPI 2009 with race (Equation 1)",
                "FIB4 units": "INDEX FIB4_ukat (AST, ALT in ukat/L as written)",
                "age": "not a covariate (paper: age is inside FIB4)",
                "missing survival outcomes": "not applied: linked mortality files are not among the downloaded public files",
                "liver exclusions": "MCQ170L == 1 (still has liver condition), MCQ510C cirrhosis or MCQ510D viral hepatitis (2017-2020 only); hepatitis serology files not downloaded",
                "PIR": "poor = PIR < 1",
                "smoking": "binary yes/no by SMQ020 (>= 100 cigarettes), Table 1 shows Smoking Yes/No",
                "drinking": "current drinker = any drinking day in past 12 months (ALQ120Q > 0; ALQ121 1-10); skipped non-drinkers = no",
                "physical activity": "low = MET-min/week < 600; GPAQ 2007+; 2005-2006 leisure only from PAQIAF",
                "education": "less than high school / high school or GED / college or above (DMDEDUC2 1-2, 3, 4-5)",
                "hypertension": "hypertension_any_140 (BP >= 140/90, self-report, medication)",
                "diabetes": "DIQ010 == 1, insulin DIQ050, pills DIQ070/DID070, FPG >= 126, random serum glucose >= 200, OGTT >= 200 (no HbA1c, as written)",
                "race": "RIDRETH1 5 groups",
                "weight": "sample weights not named; pooled MEC weight (index and outcome need the MEC exam only)"},
            "notes": "builder by agent 15/11"}

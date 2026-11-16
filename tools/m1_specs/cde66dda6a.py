"""M1 builder, cde66dda6a (METS-VF and CKD, Sci Rep 2025, doi [DOI withheld]).
Spec: METS-VF per 1 unit, CKD = eGFR < 60 or UACR > 30 (either, reading_1), Model 3, NHANES 1999-2000 to 2017-2018, weighted.
Sources in fulltext_txt/cde66dda6a.txt: population line 13, formula lines 18-19, CKD and eGFR lines 24-25, covariates lines 28-30,
weighting line 33, Table 1 (categories) line 41, Table 2 footnote line 45.
"""
import numpy as np
import pandas as pd

REC_ID = "cde66dda6a"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014",
          "2015-2016", "2017-2018"]
ALQ121_DAYS = {0: 0, 1: 365, 2: 286, 3: 182, 4: 104, 5: 52, 6: 30, 7: 12, 8: 9, 9: 4.5, 10: 1.5}
PAQIAF = {"1999-2000": "PAQIAF", "2001-2002": "PAQIAF_B", "2003-2004": "PAQIAF_C", "2005-2006": "PAQIAF_D"}


def _miss(v, codes):
    return v.where(~v.isin(codes))


def _annual_drinks(lib, df):
    """Drinks in the past 12 months = drinking days x drinks per drinking day (ALQ130); 0 when skipped as non-drinker."""
    s = lib.stems("ALQ")
    q = _miss(lib.raw(df, s, "ALQ120Q"), [777, 999])
    u = _miss(lib.raw(df, s, "ALQ120U"), [7, 9])
    days = (q * u.map({1: 52, 2: 12, 3: 1})).where(q != 0, 0.0)
    q121 = _miss(lib.raw(df, s, "ALQ121"), [77, 99])
    days = days.fillna(q121.map(ALQ121_DAYS))
    skip = (lib.raw(df, s, "ALQ110") == 2) | (lib.raw(df, s, "ALQ111") == 2)
    days = days.where(~(skip & days.isna()), 0.0)
    per = _miss(lib.raw(df, s, "ALQ130"), [77, 99, 777, 999])
    return (days * per).where(days > 0, days)


def _weekly_mv_sessions(lib, df):
    """Moderate or vigorous leisure exercise sessions per week."""
    out = pd.Series(np.nan, index=df.index)
    s = lib.stems("PAQ")
    vg = lib.raw(df, s, "PAQ650"); vd = _miss(lib.raw(df, s, "PAQ655"), [77, 99])
    mg = lib.raw(df, s, "PAQ665"); md = _miss(lib.raw(df, s, "PAQ670"), [77, 99])
    gpaq = vg.isin([1, 2]) & mg.isin([1, 2])
    sess = vd.where(vg == 1, 0.0) + md.where(mg == 1, 0.0)
    out[gpaq] = sess[gpaq]
    for c, st in PAQIAF.items():
        m = (df.cycle == c).to_numpy()
        if not m.any():
            continue
        x = lib._read(c, st)
        if x is None:
            continue
        t = _miss(x["PADTIMES"], [77777, 99999])
        wk = (t / 30 * 7).groupby(x["SEQN"]).sum()
        v = wk.reindex(df.loc[m, "SEQN"]).to_numpy(float)
        v = np.where(np.isnan(v), 0.0, v)        # PAQ respondents with no leisure activity listed
        v = np.where(df.loc[m, "has_paq"].to_numpy(bool), v, np.nan)
        out[m] = v
    return out


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["METS-VF_3.329"][0](df)
    df["y"] = ((df.egfr_ckdepi2009 < 60) | (df.uacr > 30)).astype(float).where(df.uacr.notna() & df.egfr_ckdepi2009.notna())
    df["pir3"] = pd.Series(np.select([df.pir < 1.3, df.pir <= 3.5, df.pir > 3.5], ["low", "mid", "high"], default=""),
                           index=df.index, dtype=object).replace("", np.nan)
    ad = _annual_drinks(lib, df)
    df["alcohol"] = (ad >= 12).astype(float).where(ad.notna())
    ws = _weekly_mv_sessions(lib, df)
    df["exercise"] = (ws >= 2).astype(float).where(ws.notna())
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += lib.dummies(df, "pir3", "low")
    terms += ["alcohol", "dx_hypertension", "dx_diabetes", "exercise", "albumin_serum", "scr", "ldl", "hdl", "tg",
              "u_alb", "u_cr"]
    df["w"] = lib.pooled_weight(df, "fasting")
    dom = (df.age >= 20) & df.x.notna() & df.y.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "CKD criteria": "either criterion: eGFR < 60 or UACR > 30 (reading_1)",
                "survey weight": "weighted (reading_1, variable not named); pooled fasting-subsample weight per brief rule 8 since METS-IR uses FPG and TG and the model uses LDL",
                "eGFR": "CKD-EPI 2009 with race (printed coefficients)",
                "METS-VF": "INDEX METS-VF_3.329, male = 1",
                "race": "RIDRETH1 5 groups",
                "education": "below high school / high school / above high school (Table 1)",
                "PIR": "3 categories < 1.3, 1.3-3.5, > 3.5 (Table 1)",
                "alcohol": ">= 12 drinks in past 12 months = drinking days (ALQ120Q/U or ALQ121) x ALQ130 drinks per day; skipped non-drinkers = no",
                "exercise routine": ">= 2 moderate or vigorous leisure sessions in the past week: GPAQ PAQ655 + PAQ670 days (2007+); PAQIAF times in 30 days x 7/30 (1999-2006); 'exercise plan' has no NHANES item",
                "hypertension, diabetes": "definition not given: self-reported diagnosis (BPQ020, DIQ010)",
                "continuous labs": "serum albumin g/dL, serum creatinine mg/dL, LDL, HDL, TG mg/dL, urine albumin mg/L, urine creatinine mg/dL"},
            "notes": "builder by agent 15/11"}

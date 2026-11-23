"""M1 builder, 250e5761d5 (AVAI and CKD; source withheld).
Spec: AVAI per 1 unit, CKD = eGFR < 60 or ACR > 30 (reading_1), Model 3 (Table 2 footnote list, no PIR, reading_1),
NHANES 2007-2008 to 2017-2018, unweighted (reading_1).
Formula: MathML display formula in fulltext/250e5761d5.xml (AVAI = -16.186 - 1.369 HDL + 0.38 WC + 0.144 Age - 0.013 BMI
- 0.151 TG); HDL and TG units not stated, INDEX AVAI_mmol (mmol/L, reading_1).
Sources in fulltext_txt/250e5761d5.txt: population and "No record" categories line 15, CKD line 24, covariates line 27,
Table 1 lines 120-164, Table 2 footnote line 174.
"""
import numpy as np
import pandas as pd

REC_ID = "250e5761d5"
CYCLES = ["2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]
ALQ121_DAYS = {0: 0, 1: 365, 2: 286, 3: 182, 4: 104, 5: 52, 6: 30, 7: 12, 8: 9, 9: 4.5, 10: 1.5}


def _miss(v, codes):
    return v.where(~v.isin(codes))


def _met_week(lib, df):
    """GPAQ MET-min/week: 8 x vigorous (work, leisure) + 4 x moderate (work, leisure) + 4 x walk/bike; missing parts count 0."""
    s = lib.stems("PAQ")
    tot = pd.Series(0.0, index=df.index)
    for gate, d, m, met in [("PAQ605", "PAQ610", "PAD615", 8), ("PAQ620", "PAQ625", "PAD630", 4), ("PAQ635", "PAQ640", "PAD645", 4),
                            ("PAQ650", "PAQ655", "PAD660", 8), ("PAQ665", "PAQ670", "PAD675", 4)]:
        g = lib.raw(df, s, gate)
        dd = _miss(lib.raw(df, s, d), [77, 99])
        mm = _miss(lib.raw(df, s, m), [7777, 9999])
        tot += (met * dd * mm).where(g == 1).fillna(0.0)
    return tot


def _cat(conds, labels, df, default=""):
    return pd.Series(np.select(conds, labels, default=default), index=df.index, dtype=object).replace("", np.nan)


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["AVAI_mmol"][0](df)
    df["y"] = ((df.egfr_ckdepi2009 < 60) | (df.uacr > 30)).astype(float).where(df.uacr.notna() & df.egfr_ckdepi2009.notna())
    s = lib.stems
    df["race4"] = df.race1.replace({"other_hispanic": "others", "other_multi": "others"})
    df["bmi3"] = _cat([df.bmi < 25, df.bmi <= 30, df.bmi > 30], ["normal", "over", "obese"], df)
    df["edu"] = df.educ3.astype(object)
    df["smoke"] = _cat([df.SMQ020 == 1, df.SMQ020 == 2], ["yes", "no"], df)
    # drinking: >= 12 drinks in the past year, ALQ101 (2007-2016); 2017-2018 drinking days (ALQ121) x ALQ130 >= 12
    a = lib.stems("ALQ")
    a101 = lib.raw(df, a, "ALQ101")
    q121 = _miss(lib.raw(df, a, "ALQ121"), [77, 99])
    per = _miss(lib.raw(df, a, "ALQ130"), [77, 99, 777, 999])
    yr = q121.map(ALQ121_DAYS) * per.where(q121 > 0, 0.0)
    yr = yr.where(~((lib.raw(df, a, "ALQ111") == 2) & yr.isna()), 0.0)
    df["drink"] = _cat([a101 == 1, a101 == 2, yr >= 12, yr < 12], ["yes", "no", "yes", "no"], df)
    pq = lib.stems("PAQ")
    answered = pd.concat([lib.raw(df, pq, g).isin([1, 2]) for g in ["PAQ605", "PAQ620", "PAQ635", "PAQ650", "PAQ665"]], axis=1).any(axis=1)
    df["inactive"] = (_met_week(lib, df) < 600).astype(float).where(answered)   # no GPAQ answer at all = missing (prereg 4)
    df["diab"] = _cat([df.DIQ010 == 1, df.DIQ010.isin([2, 3])], ["yes", "no"], df)
    df["htn"] = _cat([df.BPQ020 == 1, df.BPQ020 == 2], ["yes", "no"], df)
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race4", "nh_white")
    terms += lib.dummies(df, "bmi3", "normal")
    terms += lib.dummies(df, "edu", "lt_hs")
    terms += lib.dummies(df, "smoke", "no")
    terms += lib.dummies(df, "drink", "no")
    terms += ["inactive"]
    terms += lib.dummies(df, "diab", "no")
    terms += lib.dummies(df, "htn", "no")
    dom = (df.age >= 20) & df.x.notna() & df.y.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": None,
            "readings": {
                "covariate list": "Table 2 footnote list without PIR (reading_1)",
                "survey weighting": "unweighted (reading_1)",
                "ACR cut point": "> 30 mg/g (reading_1)",
                "AVAI units": "HDL and TG in mmol/L (INDEX AVAI_mmol), WC cm, BMI kg/m2, age years",
                "eGFR": "CKD-EPI 2009 with race (version not stated)",
                "age": "continuous (Table 1 shows <50/>50 only descriptively)",
                "race": "4 groups as Table 1: Mexican American, NH Black, NH White, Others (Other Hispanic + other/multi)",
                "BMI": "normal < 25, overweight 25-30 (30 inclusive), obese > 30",
                "missing covariates": "left missing, complete case (prereg section 4 overrides the paper's No record level)",
                "smoking": "SMQ020 >= 100 cigarettes yes/no",
                "drinking": "ALQ101 == 1 (>= 12 drinks/1 yr) for 2007-2016; 2017-2018 drinking days (ALQ121) x ALQ130 >= 12; ALQ111 never = no",
                "physical activity": "GPAQ MET-min/week < 600 inactive; missing domains count 0 (Table 1 has no No record level)",
                "diabetes, hypertension": "self-reported physician diagnosis DIQ010 (borderline = no), BPQ020"},
            "notes": "builder by agent 15/11. Post-run change required by prereg section 4 missing-data rule, not by n: removed the No record level of education, smoking, drinking, diabetes, hypertension; missing stays NaN and complete case drops it; physical activity is NaN when no GPAQ gate (PAQ605/620/635/650/665) was answered (main session, same rule)"}

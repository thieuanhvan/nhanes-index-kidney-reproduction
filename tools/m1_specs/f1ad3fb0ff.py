"""M1 builder, f1ad3fb0ff (Wang et al., Front Endocrinol 2023, doi [DOI withheld]). Template example.
Spec: TyG continuous, albuminuria UACR > 30 mg/g, Model 3, NHANES 2011-2018, fasting subsample weights.
Sources in the paper (fulltext_txt/f1ad3fb0ff.txt): covariates section 2.4 (lines 177-200), weights section 2.5 (line 207).
"""
import numpy as np

REC_ID = "f1ad3fb0ff"
CYCLES = ["2011-2012", "2013-2014", "2015-2016", "2017-2018"]


def build(lib):
    df = lib.table(CYCLES)
    df["x"] = lib.INDEX["TyG"][0](df)
    df["y"] = (df.uacr > 30).astype(float).where(df.uacr.notna())
    # "physical activity (vigorous/moderate/less than moderate)": GPAQ work or recreational, any vigorous > any moderate > less
    s = lambda base: lib.stems(base)
    vig = (lib.raw(df, s("PAQ"), "PAQ605") == 1) | (lib.raw(df, s("PAQ"), "PAQ650") == 1)
    mod = (lib.raw(df, s("PAQ"), "PAQ620") == 1) | (lib.raw(df, s("PAQ"), "PAQ665") == 1)
    known = lib.raw(df, s("PAQ"), "PAQ605").isin([1, 2]) | lib.raw(df, s("PAQ"), "PAQ650").isin([1, 2])
    df["pa"] = np.select([vig, mod], ["vigorous", "moderate"], default="less").astype(object)
    df.loc[~known, "pa"] = np.nan
    # "BMI was categorized as <25, 25-29.9 and >=30"
    df["bmi_cat"] = np.select([df.bmi < 25, df.bmi < 30, df.bmi >= 30], ["normal", "over", "obese"], default="").astype(object)
    df.loc[df.bmi.isna(), "bmi_cat"] = np.nan
    terms = ["x", "age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "smoking", "never")
    terms += lib.dummies(df, "pa", "less")
    terms += ["dx_hypertension", "dx_diabetes", "cvd_history"]
    terms += lib.dummies(df, "bmi_cat", "normal")
    terms += ["scr", "egfr_ckdepi2009"]
    df["w"] = lib.pooled_weight(df, "fasting")
    dom = (df.age > 20) & df.x.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {"albuminuria threshold": "UACR > 30 (reading_1)", "age": "older than 20 read as age > 20",
                         "hypertension, diabetes": "self-reported diagnosis (definition not given)",
                         "eGFR equation": "CKD-EPI 2009 with race (equation 'consisting of age, gender, race, and Scr')",
                         "race": "RIDRETH1 5 groups"},
            "notes": "template builder written by the main session"}

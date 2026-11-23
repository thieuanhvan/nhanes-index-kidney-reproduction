"""M1 builder, 259b9d6255 (WWI and CKD in hyperuricemia; source withheld).
Spec: WWI Q4 vs Q1, CKD = eGFR < 60 or albuminuria >= 30 mg/g, Model 3 of Table 2, NHANES 2007-2018.
Sources in fulltext_txt/259b9d6255.txt: population line 15, WWI line 18, hyperuricemia and CKD line 21,
covariates line 24, Model 3 line 28 and Table 2 footnote (line 145), Table 1 (lines 96-134).
"""
import numpy as np, pandas as pd

REC_ID = "259b9d6255"
CYCLES = ["2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]
URATE_LOWERING = r"ALLOPURINOL|FEBUXOSTAT|PROBENECID|LESINURAD|PEGLOTICASE|SULFINPYRAZONE"


def _flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def _bio(lib, df, col):
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for c in df.cycle.unique():
        x = lib.nc.read(c, "BIOPRO")
        if x is None or col not in x.columns:
            continue
        m = (df.cycle == c).to_numpy()
        out[m] = x.drop_duplicates("SEQN").set_index("SEQN")[col].reindex(df.loc[m, "SEQN"]).to_numpy(float)
    return out


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["WWI"][0](df).round(2)        # "rounded to two decimal places"
    egfr = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (egfr < 60)).astype(float).where(df.uacr.notna() & egfr.notna())

    df["age_g"] = np.select([df.age < 40, df.age < 65, df.age >= 65], ["young", "middle", "old"], default="").astype(object)
    df.loc[df.age.isna(), "age_g"] = np.nan
    df["smk100"] = df.SMQ020.map({1.0: 1.0, 2.0: 0.0})
    # blood pressure: hypertension = diagnosis, BP drug, mean SBP >= 130 or DBP >= 80
    df["htn"] = _flag((df.BPQ020 == 1) | (df.BPQ_med_now == 1) | (df.sbp >= 130) | (df.dbp >= 80),
                      df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna()))
    # blood glucose: diabetes = diagnosis, insulin or pills, HbA1c > 6.5, FBG > 126 (strict as written)
    pills = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    pos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (pills == 1) | (df.hba1c > 6.5) | (df.glucose > 126)
    df["dm"] = _flag(pos, df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna())
    # blood lipid: hyperlipidemia = diagnosis, lipid-lowering drug, TC >= 240, HDL < 40/50, LDL >= 160, TG >= 150
    told = lib.raw(df, s("BPQ"), "BPQ080"); med = lib.raw(df, s("BPQ"), "BPQ100D")
    hdl_low = ((df.female == 0) & (df.hdl < 40)) | ((df.female == 1) & (df.hdl < 50))
    pos = (told == 1) | (med == 1) | (df.tc >= 240) | hdl_low | (df.ldl >= 160) | (df.tg >= 150)
    df["hyperlip"] = _flag(pos, told.isin([1, 2]) | df.tc.notna() | df.hdl.notna())
    df["phos"] = _bio(lib, df, "LBXSPH")
    df["ult"] = lib.rx_any(df, URATE_LOWERING)

    ua = df.uric_acid
    hyperua = ((df.female == 0) & (ua > 7)) | ((df.female == 1) & (ua > 6))
    dom = (hyperua & (df.age >= 20) & df.waist.notna() & df.weight.notna() & df.ult.notna() & df.x.notna()
           & df.uacr.notna() & egfr.notna())
    qn, cuts = lib.quantile_dummies(df, "x", dom, 4, prefix="x_q")
    terms = qn
    terms += lib.dummies(df, "age_g", "young")
    terms += ["female"]
    terms += lib.dummies(df, "race1", "mexican_american")
    terms += lib.dummies(df, "educ3", "lt_hs")
    terms += ["pir", "smk100", "htn", "dm", "hyperlip", "calcium", "phos"]

    df["w"] = lib.pooled_weight(df, "mec")
    readings = {
        "exposure scale": "Q4 vs Q1 (reading_1); quartiles cut unweighted on the domain (rule 6); paper's printed cuts 10.72, 11.27, 11.84",
        "UACR threshold": "albuminuria >= 30 mg/g (reading_1)",
        "eGFR": "CKD-EPI 2009 with race (version not stated)",
        "hyperuricemia": "serum uric acid above 7 mg/dL men, above 6 mg/dL women (strict, 'above')",
        "urate-lowering medication": "missing = no RXQ_RX record; drug use by name (allopurinol, febuxostat, probenecid, lesinurad, pegloticase, sulfinpyrazone); used only for the exclusion as the paper lists it",
        "WWI": "waist (cm) / sqrt(weight kg), rounded to 2 decimals as written",
        "blood pressure, blood glucose, blood lipid": "the three binary conditions of Methods (4) and Table 1 (Pressure, Glucose, Lipids)",
        "hypertension": "BPQ020=1 or BPQ050A=1 or mean SBP >= 130 or DBP >= 80",
        "diabetes": "DIQ010=1, insulin DIQ050=1, pills DIQ070/DID070, HbA1c > 6.5, fasting glucose > 126 (strict as written)",
        "hyperlipidemia": "BPQ080=1, BPQ100D=1, TC >= 240, HDL < 40 men / < 50 women, LDL >= 160, TG >= 150 (fasting TRIGLY, criterion evaluated when present)",
        "age": "3 groups 20-39 / 40-64 / >= 65 as Methods define",
        "PIR": "continuous: Methods name PIR without categories; Table 1 groups have no stated cut points",
        "smoking": ">= 100 cigarettes yes/no",
        "education, race": "3 levels; RIDRETH1 5 groups",
        "calcium, phosphorus": "serum, mg/dL, continuous",
        "survey weighting": "weighted, MEC weight (no covariate requires the fasting subsample)",
    }
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x_q4", "weight": "w",
            "readings": readings, "notes": "builder written by an M1 agent 15/11/2026; quartile cuts on the domain " +
            ", ".join(str(int(c * 100) / 100) for c in cuts)}

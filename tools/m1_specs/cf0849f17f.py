"""M1 builder, cf0849f17f (SII and DKD in T2DM, Front Endocrinol 2022, doi [DOI withheld]).
Spec: SII >= 445.21 vs < 445.21 (fixed cut point), DKD = ACR >= 30 mg/g and/or eGFR < 60 in T2DM, Model 3 of Table 2,
NHANES 2011-2018. Sources in fulltext_txt/cf0849f17f.txt: population lines 246-252, diabetes and DKD lines 260-267,
other definitions lines 269-297, Model 3 lines 336-341 and Table 2 footnote (line ~720).
"""
import numpy as np, pandas as pd

REC_ID = "cf0849f17f"
CYCLES = ["2011-2012", "2013-2014", "2015-2016", "2017-2018"]

NSAID_STEROID = (r"IBUPROFEN|NAPROXEN|CELECOXIB|DICLOFENAC|MELOXICAM|INDOMETHACIN|KETOROLAC|ETODOLAC|NABUMETONE|"
                 r"PIROXICAM|SULINDAC|OXAPROZIN|KETOPROFEN|FLURBIPROFEN|DIFLUNISAL|SALSALATE|MEFENAMIC|FENOPROFEN|"
                 r"TOLMETIN|ASPIRIN|PREDNISONE|PREDNISOLONE|METHYLPREDNISOLONE|DEXAMETHASONE|HYDROCORTISONE|CORTISONE")


def _flag(pos, obs):
    return pd.Series(np.where(pos, 1.0, np.where(obs, 0.0, np.nan)), index=pos.index)


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    sii = lib.INDEX["SII"][0](df)
    df["sii"] = sii
    df["x"] = (sii >= 445.21).astype(float).where(sii.notna())
    egfr = df.egfr_ckdepi2009
    df["y"] = ((df.uacr >= 30) | (egfr < 60)).astype(float).where(df.uacr.notna() & egfr.notna())

    # diabetes (lines 260-263): diagnosis, FPG >= 7.0 mmol/L, HbA1c >= 6.5, diabetes drugs
    pos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (df.DIQ070 == 1) | (df.hba1c >= 6.5) | (df.glucose >= 126)
    obs = df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna()
    df["dm"] = _flag(pos, obs)
    # hypertension (lines 269-271): SBP >= 140 and/or DBP >= 90 or physician diagnosis
    df["htn"] = _flag((df.sbp >= 140) | (df.dbp >= 90) | (df.BPQ020 == 1), df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna()))
    # hyperlipidemia (lines 271-274): any of TC >= 240, TG >= 200, LDL >= 160, HDL < 40, or diagnosis (BPQ080)
    told = lib.raw(df, s("BPQ"), "BPQ080")
    pos = (df.tc >= 240) | (df.tg_serum >= 200) | (df.ldl >= 160) | (df.hdl < 40) | (told == 1)
    df["hyperlip"] = _flag(pos, df.tc.notna() | df.hdl.notna() | told.isin([1, 2]))
    # MetS, NCEP ATP III (revised): >= 3 of waist, TG, HDL, BP, fasting glucose; diabetics not excluded
    bpmed = df.BPQ_med_now == 1
    dmtx = (df.DIQ050 == 1) | (df.DIQ070 == 1)
    comp = pd.DataFrame({
        "waist": np.where(df.waist.isna(), np.nan, ((df.female == 1) & (df.waist > 88)) | ((df.female == 0) & (df.waist > 102))),
        "tg": np.where(df.tg_serum.isna(), np.nan, df.tg_serum >= 150),
        "hdl": np.where(df.hdl.isna(), np.nan, ((df.female == 1) & (df.hdl < 50)) | ((df.female == 0) & (df.hdl < 40))),
        "bp": np.where(bpmed | (df.sbp >= 130) | (df.dbp >= 85), 1.0, np.where(df.sbp.notna() & df.dbp.notna(), 0.0, np.nan)),
        "glu": np.where(dmtx | (df.glucose >= 100), 1.0, np.where(df.glucose.notna(), 0.0, np.nan)),
    }, index=df.index).astype(float)
    npos = comp.sum(axis=1); nmiss = comp.isna().sum(axis=1)
    df["mets"] = np.where(npos >= 3, 1.0, np.where(npos + nmiss < 3, 0.0, np.nan))
    df["antiinfl"] = lib.rx_any(df, NSAID_STEROID)
    df["metformin"] = lib.rx_any(df, r"METFORMIN")
    df["insulin_use"] = df.DIQ050.map({1.0: 1.0, 2.0: 0.0})
    mcq = lambda c: lib.raw(df, s("MCQ"), c)
    m = {c: mcq(c) for c in ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]}
    asc = ["MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]
    df["ascvd"] = _flag(pd.concat([m[c] == 1 for c in asc], axis=1).any(axis=1),
                        pd.concat([m[c].isin([1, 2]) for c in asc], axis=1).all(axis=1))
    df["chf"] = m["MCQ160B"].map({1.0: 1.0, 2.0: 0.0})
    df["bmi_g"] = np.select([df.bmi < 25, df.bmi < 30, df.bmi >= 30], ["normal", "over", "obese"], default="").astype(object)
    df.loc[df.bmi.isna(), "bmi_g"] = np.nan

    terms = ["x", "age", "female", "pir"]
    terms += lib.dummies(df, "race1", "mexican_american")
    terms += lib.dummies(df, "educ3", "gt_hs")
    terms += lib.dummies(df, "smoking", "former")
    terms += ["htn", "hyperlip", "mets", "antiinfl", "metformin", "insulin_use"]
    terms += lib.dummies(df, "bmi_g", "normal")
    terms += ["sbp", "dbp", "hba1c", "scr", "bun", "uric_acid", "egfr_ckdepi2009", "tc", "tg_serum", "alt", "ast",
              "ascvd", "chf"]

    df["w"] = lib.pooled_weight(df, "mec")
    dom = (df.age >= 18) & (df.pregnant != 1) & (df.dm == 1) & df.sii.notna() & df.uacr.notna() & egfr.notna()
    readings = {
        "Model 3 covariates": "Table 2 footnote list incl. BMI, serum creatinine and eGFR (reading_1)",
        "T2DM": "diabetes by the paper's four criteria (no type 1 exclusion is stated); borderline not diabetes; FPG >= 126 mg/dL, HbA1c >= 6.5 %",
        "eGFR": "CKD-EPI 2009 with race (version not stated)",
        "UACR threshold": ">= 30 mg/g as written",
        "exposure": "SII >= 445.21 vs < 445.21, counts in 10^9/L as in NHANES",
        "survey weighting": "weighted, MEC weight (CBC, biochemistry, urine; no fasting restriction in the flow)",
        "hypertension": "SBP >= 140 or DBP >= 90 (mean of readings) or BPQ020 = 1; medication not part of the definition",
        "hyperlipidemia": "any one criterion (the 'and' before HDL read as a list), physician diagnosis = BPQ080 'told high cholesterol'; LDL fasting subsample only, criterion evaluated when present",
        "TG": "serum biochemistry TG (LBXSTR), the paper lists TG in the biochemical profile with Scr, BUN, SUA, ALT, AST",
        "MetS": "NCEP ATP III revised: waist > 102/88 cm, TG >= 150, HDL < 40/50, BP >= 130/85 or BP drug, FPG >= 100 or insulin/pills; MetS = >= 3 positive, 0 when fewer than 3 could be positive",
        "anti-inflammation therapy": "any prescription (RXQ_RX) of an NSAID, aspirin or systemic glucocorticoid by drug name",
        "metformin use": "any prescription containing metformin (RXQ_RX)",
        "insulin use": "DIQ050 taking insulin now",
        "ASCVD": "self-reported CHD, angina, MI or stroke; CHF: MCQ160B",
        "covariate coding": "age, PIR, SBP, DBP, labs, eGFR continuous; BMI WHO 3 groups (paper's own classification); smoking never/former/now; education 3 levels, 'Others' kept missing; race RIDRETH1 5 groups",
    }
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": readings, "notes": "builder written by an M1 agent 15/11/2026"}

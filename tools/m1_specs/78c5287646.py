"""M1 builder, 78c5287646 (PNI and DKD in T2DM; source withheld).
Spec: PNI > 50.5 vs <= 50.5, DKD = UACR > 30 mg/g and/or eGFR < 60 in T2DM, multivariate model of Table 5,
NHANES 2013-2018, 6-year weight = 2-year weight / 3. Sources in fulltext_txt/78c5287646.txt: population line 12,
PNI line 16, definitions lines 19-20, weights line 23, Table 5 and footnote lines 181-202.
"""
import numpy as np, pandas as pd
import diet_indices_i139 as di

REC_ID = "78c5287646"
CYCLES = ["2013-2014", "2015-2016", "2017-2018"]


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
    pni = lib.INDEX["PNI"][0](df)            # 5 x lymphocytes (10^9/L) + albumin (g/L)
    df["pni"] = pni
    df["x"] = (pni > 50.5).astype(float).where(pni.notna())
    egfr = df.egfr_ckdepi2009
    df["y"] = ((df.uacr > 30) | (egfr < 60)).astype(float).where(df.uacr.notna() & egfr.notna())

    # diabetes (line 19): self-report, medication or insulin, HbA1c >= 6.5, FPG >= 7.0, random glucose >= 11.1 mmol/L
    glu_si = _bio(lib, df, "LBDSGLSI")
    pos = ((df.DIQ010 == 1) | (df.DIQ050 == 1) | (df.DIQ070 == 1) | (df.hba1c >= 6.5) | (df.glucose >= 126)
           | (glu_si >= 11.1))
    obs = df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna() | glu_si.notna()
    df["dm"] = _flag(pos, obs)
    # type 1 (line 12): diagnosed before age 30 and insulin as the only hypoglycemic agent
    age_dx = lib.raw(df, s("DIQ"), "DID040").where(lambda v: v < 100)
    t1 = (age_dx < 30) & (df.DIQ050 == 1) & (df.DIQ070 != 1)

    # covariates of the Table 5 footnote, coded as Table 5 shows them
    df["male"] = 1 - df.female
    df["nhw"] = (df.race1 == "nh_white").astype(float).where(df.race1.notna())
    df["ever_smk"] = df.smoking.map({"never": 0.0, "former": 1.0, "current": 1.0})
    df["htn"] = _flag((df.BPQ020 == 1) | (df.BPQ_med_now == 1) | (df.sbp > 140) | (df.dbp > 90),
                      df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna()))
    d1 = di.diet(df, "day1")
    dii = di.dii(d1, di.DII28)
    df["dii_pos"] = (dii > 0).astype(float).where(dii.notna())
    wt = df.weight.where(df.weight > 0)
    df["energy_kg"] = d1["kcal"] / wt
    df["dpi"] = d1["prot"] / wt
    ua = _bio(lib, df, "LBDSUASI")
    df["hyperua"] = _flag(((df.female == 0) & (ua >= 420)) | ((df.female == 1) & (ua >= 360)), ua.notna() & df.female.notna())
    df["tg_mmol"] = _bio(lib, df, "LBDSTRSI")
    alb = df.albumin_serum * 10
    df["hypoalb"] = (alb < 35).astype(float).where(alb.notna())
    hgb = df.hgb * 10
    df["anemia"] = _flag(((df.female == 0) & (hgb < 130)) | ((df.female == 1) & (hgb < 120)), hgb.notna() & df.female.notna())

    terms = ["x", "age", "male", "nhw", "ever_smk", "htn", "dii_pos", "energy_kg", "dpi", "hyperua", "hba1c", "tg_mmol",
             "hypoalb", "anemia", "lymph_n"]

    w1 = lib.raw(df, s("DR1TOT"), "WTDRD1")
    df["w"] = (w1 / 3).fillna(0)
    dom = ((df.age >= 18) & (df.pregnant != 1) & (df.dm == 1) & ~t1 & d1["kcal"].notna() & df.hba1c.notna()
           & egfr.notna() & df.uacr.notna() & df.lymph_n.notna() & df.pni.notna())
    readings = {
        "analytic N": "the T2DM participants (reading_1)",
        "covariates": "footnote list without the exposure 'PNI group' (reading_1)",
        "inflection point": "50.5 used as the fixed cut point, PNI > 50.5 vs <= 50.5",
        "eGFR": "CKD-EPI 2009 with race (CKD_EPI_Scr, version not stated)",
        "UACR threshold": "> 30 mg/g as written",
        "diabetes": "DIQ010=1, insulin DIQ050=1, pills DIQ070=1, HbA1c >= 6.5, fasting glucose >= 126 mg/dL, serum glucose (LBDSGLSI) >= 11.1 mmol/L as random glucose; borderline not diabetes",
        "type 1 exclusion": "age at diagnosis DID040 < 30 and DIQ050 = 1 and not DIQ070 = 1",
        "mortality data exclusion": "not applied: the NCHS linked mortality file is not in the downloaded CDC NHANES files; only adults with too little identifying data are ineligible",
        "dietary data": "valid day-1 recall (DR1DRSTZ = 1) required, the first-day recall is the only one used",
        "DII": "Shivappa protocol with the usual 28 NHANES parameters (the paper names no parameter list), day-1 recall, constants as in tools/diet_indices_i139.py; DII group = DII > 0 vs <= 0 (Table 1 and 5 groups 'DII < 0' / 'DII > 0')",
        "energy intake": "kcal/kg/d (Table 5 unit, the model table): day-1 energy / measured body weight",
        "DPI": "day-1 protein (g) / measured body weight (kg)",
        "race": "Non-Hispanic White yes vs no (Table 5 coding)",
        "smoking": "ever (former or current) vs never (Table 5 'Smoking (yes vs. no)')",
        "hypertension": "self-report BPQ020=1, BPQ050A=1, mean SBP > 140 or DBP > 90 (strict, as written)",
        "hyperuricemia": "LBDSUASI >= 420 umol/L men, >= 360 women",
        "triglyceride": "serum biochemistry TG in mmol/L (LBDSTRSI): Table 1 reports TG for all 14,349 and the exclusion list has no fasting restriction",
        "hypoalbuminemia, anemia": "albumin < 35 g/L; hemoglobin < 130 g/L men, < 120 women",
        "survey weighting": "dietary day-1 weight WTDRD1 / 3 (the paper uses dietary data, the smallest subsample; one third of the 2-year weight as stated)",
    }
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": readings, "notes": "builder written by an M1 agent 15/11/2026"}

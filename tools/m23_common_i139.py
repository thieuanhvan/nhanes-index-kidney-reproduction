"""I139 M2/M3 shared code, following extraction/m2_m3_rules_i139.md (locked before any M2/M3 run).
Loads a builder (tools/m1_specs/<rec_id>.py), the correct weight, standardised exposure, concept codings for added
covariates, the outcome-component flag, Rao-Wu replicate weights and a fast weighted IRLS point fit.
"""
import os, sys, re, json, hashlib, importlib.util
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m1_lib_i139 as lib
import nhanes_core as nc
import nhanes_extra_i139 as nx

ROOT = os.path.dirname(HERE)
X = os.path.join(ROOT, "extraction")
LOG110 = float(np.log(1.10))
ALQ121_DAYS = {0: 0, 1: 365, 2: 286, 3: 182, 4: 104, 5: 52, 6: 30, 7: 12, 8: 9, 9: 4.5, 10: 1.5}


def spec_table():
    lock = open(os.path.join(X, "spec_table_i139.sha256"), encoding="utf-8").read()
    pj = os.path.join(X, "spec_table_i139.json")
    assert hashlib.sha256(open(pj, "rb").read()).hexdigest() in lock, "spec table changed after lock"
    return {r["rec_id"]: r for r in json.load(open(pj, encoding="utf-8"))}


def labels():
    import csv
    return {r["rec_id"]: r for r in csv.DictReader(open(os.path.join(X, "labels_i139.tsv"), encoding="utf-8"), delimiter="\t")}


def load_builder(rid):
    p = os.path.join(HERE, "m1_specs", rid + ".py")
    sp = importlib.util.spec_from_file_location("b_" + rid, p)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m.build(lib)


# ------------------------------------------------------------------ flag (rule 3)
FLAG_RX = re.compile(r"\balbumin\b|serum albumin|\bALB\b|hemoglobin|haemoglobin|\bHGB\b|\bHb\b|creatinine|cystatin|\bBUN\b|"
                     r"urea nitrogen|phosph|potassium", re.I)
UA_RX = re.compile(r"uric|\bSUA\b|\bUA\b", re.I)


def outcome_flag(spec):
    comps = " | ".join(spec["index_components"] or [])
    comps = re.sub(r"urinary albumin", "URINALB", comps, flags=re.I)
    hit = bool(FLAG_RX.search(comps)) or "URINALB" in comps
    out = (spec["outcome_main"] or "").lower()
    has_egfr = any(k in out for k in ["ckd", "dkd", "egfr", "nephropathy"]) or out.startswith("dn")
    if UA_RX.search(comps) and has_egfr:
        hit = True
    return hit


# ------------------------------------------------------------------ weights (rule 1)
def diet_weight(df):
    cyc = list(dict.fromkeys(df.cycle))
    four = "1999-2000" in cyc and "2001-2002" in cyc
    total = sum(nc.CYCLES[c][4] for c in cyc)
    w = pd.Series(0.0, index=df.index)
    st = nx.STEMS["DIET1"]
    w1 = lib.raw(df, st, ["WTDRD1", "WTDRD1PP"])
    w4 = lib.raw(df, st, "WTDR4YR") if four else None
    for c in cyc:
        m = df.cycle == c
        if c == "2017-2020":
            w[m] = w1[m].fillna(0) * 3.2 / total
        elif four and c in ("1999-2000", "2001-2002"):
            w[m] = w4[m].fillna(0) * 4 / total
        else:
            w[m] = w1[m].fillna(0) * 2 / total
    return w


def correct_weight(rid, b):
    """Return (weight Series, kind) by rule 1."""
    df = b["df"]
    if rid in lib.DIET:
        return diet_weight(df), "diet_day1"
    dom = b["domain"].fillna(False).astype(bool)
    cc = dom & df[[b["y"]] + b["terms"]].notna().all(axis=1)
    wf = lib.pooled_weight(df, "fasting")
    if cc.sum() and (wf[cc] > 0).mean() >= float(os.environ.get("I139_W_THRESHOLD", "0.80")):   # rule 1 as amended 15/11 (was 0.99)
        return wf, "fasting"
    return lib.pooled_weight(df, "mec"), "mec"


def wrong_weight(df, kind):
    return lib.pooled_weight(df, "interview" if kind == "mec" else "mec")


# ------------------------------------------------------------------ concept codings (rule 4)
def _miss(v, codes):
    return v.where(~v.isin(codes))


def add_concept(df, concept):
    """Add the fixed coding of a concept to df; return the list of term names."""
    s = lib.stems
    if concept == "age":
        return ["age"]
    if concept == "sex":
        return ["female"]
    if concept == "race":
        return lib.dummies(df, "race1", "nh_white")
    if concept == "education":
        return lib.dummies(df, "educ3", "lt_hs")
    if concept == "pir":
        return ["pir"]
    if concept == "smoking":
        return lib.dummies(df, "smoking", "never")
    if concept == "hypertension_dx":
        df["c_htn_dx"] = df.BPQ020.map({1: 1.0, 2: 0.0})
        return ["c_htn_dx"]
    if concept == "diabetes_dx":
        pills = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
        pos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (pills == 1)
        neg = df.DIQ010.isin([2, 3])
        df["c_dm_dx"] = pd.Series(np.where(pos, 1.0, np.where(neg, 0.0, np.nan)), index=df.index)
        return ["c_dm_dx"]
    if concept == "physical_activity":
        p = s("PAQ")
        g650, g665 = lib.raw(df, p, "PAQ650"), lib.raw(df, p, "PAQ665")
        p200, p320 = lib.raw(df, p, "PAD200"), lib.raw(df, p, "PAD320")
        new = df.cycle >= "2007-2008"
        a, b = g650.where(new, p200), g665.where(new, p320)
        pos = (a == 1) | (b == 1)
        neg = a.isin([2, 3]) & b.isin([2, 3])   # PAD200/PAD320 3 = unable to do activity
        df["c_pa"] = pd.Series(np.where(pos, 1.0, np.where(neg, 0.0, np.nan)), index=df.index)
        return ["c_pa"]
    if concept == "alcohol":
        a = s("ALQ")
        q = _miss(lib.raw(df, a, "ALQ120Q"), [777, 999])
        q121 = _miss(lib.raw(df, a, "ALQ121"), [77, 99])
        any_ = (q > 0).astype(float).where(q.notna())
        any_ = any_.fillna((q121 > 0).astype(float).where(q121.notna()))
        skip = (lib.raw(df, a, "ALQ110") == 2) | (lib.raw(df, a, "ALQ111") == 2)
        df["c_alc"] = any_.where(~(skip & any_.isna()), 0.0)
        return ["c_alc"]
    if concept == "lipid_drugs":
        bp = s("BPQ")
        d = lib.raw(df, bp, "BPQ100D")
        neg = (d == 2) | (lib.raw(df, bp, "BPQ090D") == 2) | (lib.raw(df, bp, "BPQ080") == 2) | (lib.raw(df, bp, "BPQ060") == 2)
        df["c_lipid_rx"] = pd.Series(np.where(d == 1, 1.0, np.where(neg, 0.0, np.nan)), index=df.index)
        return ["c_lipid_rx"]
    if concept == "bmi_waist":
        return ["bmi"]
    raise KeyError(concept)


def age_terms(terms):
    return [t for t in terms if t == "age" or t.startswith("age") and not t.startswith("age_at")]


# ------------------------------------------------------------------ standardise (rule 2)
def standardise(df, col, mask, w):
    v = df.loc[mask, col].to_numpy(float)
    if np.unique(v[~np.isnan(v)]).size <= 2:
        return col, 1.0
    ww = w[mask].to_numpy(float)
    m = np.average(v, weights=ww)
    sd = float(np.sqrt(np.average((v - m) ** 2, weights=ww)))
    df[col + "__sd"] = (df[col] - m) / sd
    return col + "__sd", sd


# ------------------------------------------------------------------ fit + Rao-Wu (rule 7)
def irls(Xm, y, w, beta0=None, it=30, tol=1e-9):
    beta = np.zeros(Xm.shape[1]) if beta0 is None else beta0.copy()
    w = w / w[w > 0].mean()
    for _ in range(it):
        eta = np.clip(Xm @ beta, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        H = Xm.T @ (Xm * (w * p * (1 - p))[:, None])
        g = Xm.T @ (w * (y - p))
        step = np.linalg.solve(H, g)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    return beta, H


def raowu_multipliers(df_design, B, seed):
    """df_design: DataFrame with SDMVSTRA, SDMVPSU of all design rows. Returns (keys array, B x K multiplier matrix)."""
    des = df_design[["SDMVSTRA", "SDMVPSU"]].drop_duplicates().sort_values(["SDMVSTRA", "SDMVPSU"])
    keys = des.SDMVSTRA.to_numpy(np.int64) * 1000 + des.SDMVPSU.to_numpy(np.int64)
    hs = des.SDMVSTRA.to_numpy()
    rng = np.random.default_rng(int(hashlib.sha256(seed.encode()).hexdigest()[:15], 16))
    Mx = np.zeros((B, len(keys)))
    for h in np.unique(hs):
        idx = np.where(hs == h)[0]
        n = len(idx)
        if n < 2:
            raise ValueError(f"stratum {h} has {n} PSU")
        for bb in range(B):
            draw = rng.integers(0, n, n - 1)
            cnt = np.bincount(draw, minlength=n)
            Mx[bb, idx] = cnt * n / (n - 1)
    return keys, Mx

"""nhanes_core: harmonised person-level NHANES table, 1999-2000 to 2021-2023, for study I139 (prereg_i139_v7_en.md, sections 4 and 6).

Infrastructure only. This module builds variables and survey-design columns and provides the Taylor-linearisation
survey logistic estimator used by the study; it does not choose any index x kidney model.

Sources: public CDC/NCHS XPT files in <workdir>/data-external/nhanes-extra/<cycle>/
(sha256 in download_log.json). Variable sources per cycle: tools/nhanes_core_variables.md (regenerate with --doc).

Rules implemented
- Missing codes 7/9 (77/99, 777/999 where the codebook uses them) of questionnaire items stay missing; no "missing" level.
- Serum creatinine calibrated to the standard (NHANES recommendation, LAB18 and BIOPRO_D documentation):
    1999-2000: Scr_std = 1.013 * Scr + 0.147 ; 2005-2006: Scr_std = -0.016 + 0.978 * Scr. Other cycles unchanged.
- eGFR CKD-EPI 2009 (race term: non-Hispanic Black, RIDRETH1 == 4) and CKD-EPI 2021 (race-free).
- UACR (mg/g) = urine albumin (mg/L) / urine creatinine (mg/dL) * 100.
- Blood pressure: mean of available readings; 2017-2018, 2017-March 2020 and 2021-2023 use the oscillometric
  protocol (BPXO*), earlier cycles the auscultatory one (BPX*); diastolic 0 is treated as missing.
- Pooled weights: weight_i * years_i / sum(years) over the pooled releases; years = 2 for a 2-year release,
  3.2 for 2017-March 2020 (WT*PRP), and the 4-year weights (WT*4YR) are used for 1999-2002 when both 1999-2000
  and 2001-2002 are pooled (years = 4 for that pair). 2017-2018 (J) and 2017-March 2020 (P) overlap and cannot be pooled.

Usage
    import nhanes_core as nc
    df = nc.build(["2005-2006", "2007-2008"])          # person-level table, all ages
    df["w"] = nc.pooled_weight(df, "fasting")          # or "mec", "interview"
    res = nc.svy_logit(df, y="...", x=[...], w="w", domain=mask)
CLI: python -I nhanes_core.py --doc     (rewrites nhanes_core_variables.md from the XPT column catalogue)
     python -I nhanes_core.py --summary (counts per cycle, no outcome modelling)
"""
import ctypes, os, sys
import numpy as np, pandas as pd
from scipy import stats

try: ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)  # below normal
except Exception: pass

ROOT = os.environ.get("NHANES_ROOT", "<workdir>/data-external/nhanes-extra")
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- cycle registry
# component -> file stem; prefix/suffix rules differ in 1999-2004 and in the pre-pandemic file
_STD = {"DEMO": "DEMO", "ALB_CR": "ALB_CR", "BIOPRO": "BIOPRO", "GLU": "GLU", "TRIGLY": "TRIGLY", "HDL": "HDL", "TCHOL": "TCHOL",
        "GHB": "GHB", "INS": "INS", "BMX": "BMX", "BPX": "BPX", "BPXO": "BPXO", "DIQ": "DIQ", "BPQ": "BPQ", "SMQ": "SMQ", "MCQ": "MCQ", "CBC": "CBC", "FASTQX": "FASTQX"}
_OLD = {"DEMO": "DEMO", "ALB_CR": "L16", "BIOPRO": "L40", "GLU": "L10AM", "TRIGLY": "L13AM", "HDL": "L13", "TCHOL": "L13",
        "GHB": "L10", "BMX": "BMX", "BPX": "BPX", "DIQ": "DIQ", "BPQ": "BPQ", "SMQ": "SMQ", "MCQ": "MCQ", "CBC": "L25"}
_9900 = {"DEMO": "DEMO", "ALB_CR": "LAB16", "BIOPRO": "LAB18", "GLU": "LAB10AM", "TRIGLY": "LAB13AM", "HDL": "LAB13", "TCHOL": "LAB13",
         "GHB": "LAB10", "BMX": "BMX", "BPX": "BPX", "DIQ": "DIQ", "BPQ": "BPQ", "SMQ": "SMQ", "MCQ": "MCQ", "CBC": "LAB25"}
CYCLES = {  # name: (start year, file map, prefix, suffix, years, oscillometric BP)
    "1999-2000": (1999, _9900, "", "", 2, False),
    "2001-2002": (2001, _OLD, "", "_B", 2, False),
    "2003-2004": (2003, _OLD, "", "_C", 2, False),
    "2005-2006": (2005, _STD, "", "_D", 2, False),
    "2007-2008": (2007, _STD, "", "_E", 2, False),
    "2009-2010": (2009, _STD, "", "_F", 2, False),
    "2011-2012": (2011, _STD, "", "_G", 2, False),
    "2013-2014": (2013, _STD, "", "_H", 2, False),
    "2015-2016": (2015, _STD, "", "_I", 2, False),
    "2017-2018": (2017, _STD, "", "_J", 2, True),
    "2017-2020": (2017, _STD, "P_", "", 3.2, True),
    "2021-2023": (2021, _STD, "", "_L", 2, True),
}
ALL_CYCLES = list(CYCLES)

def file_path(cycle, comp):
    _, fmap, pre, suf, _, _ = CYCLES[cycle]
    if comp not in fmap: return None
    p = os.path.join(ROOT, cycle, f"{pre}{fmap[comp]}{suf}.xpt")
    return p if os.path.exists(p) else None

_cache = {}
def read(cycle, comp):
    p = file_path(cycle, comp)
    if p is None: return None
    if p not in _cache:
        x = pd.read_sas(p, format="xport")
        x.columns = [c.upper() for c in x.columns]
        # pandas reads an XPT zero stored as the smallest IBM float (5.4e-79); set it to an exact 0
        num = x.select_dtypes("number").columns
        x[num] = x[num].mask(x[num].abs() < 1e-70, 0.0)
        x["SEQN"] = x["SEQN"].astype(np.int64)
        _cache[p] = x
    return _cache[p]

# ---------------------------------------------------------------- variable catalogue
# harmonised name -> (component, candidate source columns in priority order, missing codes)
QCODES = (7, 9, 77, 99, 777, 999, 7777, 9999)
VARS = {
    # demographics and design
    "age": ("DEMO", ["RIDAGEYR"], ()), "sex": ("DEMO", ["RIAGENDR"], ()), "RIDRETH1": ("DEMO", ["RIDRETH1"], ()),
    "RIDRETH3": ("DEMO", ["RIDRETH3"], ()), "educ_raw": ("DEMO", ["DMDEDUC2"], (7, 9)), "pir": ("DEMO", ["INDFMPIR"], ()),
    "pregnant_raw": ("DEMO", ["RIDEXPRG"], ()), "SDMVSTRA": ("DEMO", ["SDMVSTRA"], ()), "SDMVPSU": ("DEMO", ["SDMVPSU"], ()),
    "WTINT2YR": ("DEMO", ["WTINT2YR"], ()), "WTMEC2YR": ("DEMO", ["WTMEC2YR"], ()), "WTINT4YR": ("DEMO", ["WTINT4YR"], ()),
    "WTMEC4YR": ("DEMO", ["WTMEC4YR"], ()), "WTINTPRP": ("DEMO", ["WTINTPRP"], ()), "WTMECPRP": ("DEMO", ["WTMECPRP"], ()),
    # questionnaire
    "SMQ020": ("SMQ", ["SMQ020"], QCODES), "SMQ040": ("SMQ", ["SMQ040"], QCODES),
    "DIQ010": ("DIQ", ["DIQ010"], QCODES), "DIQ050": ("DIQ", ["DIQ050"], QCODES), "DIQ070": ("DIQ", ["DIQ070"], QCODES),
    "BPQ020": ("BPQ", ["BPQ020"], QCODES), "BPQ_med_now": ("BPQ", ["BPQ050A", "BPQ150"], QCODES),
    "MCQ160B": ("MCQ", ["MCQ160B", "MCQ160b"], QCODES), "MCQ160C": ("MCQ", ["MCQ160C", "MCQ160c"], QCODES),
    "MCQ160D": ("MCQ", ["MCQ160D", "MCQ160d"], QCODES), "MCQ160E": ("MCQ", ["MCQ160E", "MCQ160e"], QCODES),
    "MCQ160F": ("MCQ", ["MCQ160F", "MCQ160f"], QCODES),
    # examination
    "bmi": ("BMX", ["BMXBMI"], ()), "waist": ("BMX", ["BMXWAIST"], ()),
    # laboratory
    "glucose": ("GLU", ["LBXGLU"], ()), "WTSAF2YR": ("GLU", ["WTSAF2YR"], ()), "WTSAF4YR": ("GLU", ["WTSAF4YR"], ()),
    "WTSAFPRP": ("GLU", ["WTSAFPRP"], ()), "fast_hours": ("GLU", ["PHAFSTHR"], ()), "fast_hours_fq": ("FASTQX", ["PHAFSTHR"], ()),
    "WTPH2YR": ("BIOPRO", ["WTPH2YR"], ()),
    "insulin": ("INS", ["LBXIN"], ()), "insulin_glu": ("GLU", ["LBXIN"], ()),
    "hba1c": ("GHB", ["LBXGH"], ()),
    "tg": ("TRIGLY", ["LBXTR", "LBXTLG"], ()), "ldl": ("TRIGLY", ["LBDLDL"], ()),
    "hdl": ("HDL", ["LBDHDD", "LBXHDD", "LBDHDL"], ()), "tc": ("TCHOL", ["LBXTC"], ()),
    "scr_raw": ("BIOPRO", ["LBXSCR", "LBDSCR"], ()), "tg_serum": ("BIOPRO", ["LBXSTR"], ()), "glucose_serum": ("BIOPRO", ["LBXSGL"], ()),
    "tc_serum": ("BIOPRO", ["LBXSCH"], ()), "uric_acid": ("BIOPRO", ["LBXSUA"], ()), "alt": ("BIOPRO", ["LBXSATSI"], ()),
    "ast": ("BIOPRO", ["LBXSASSI"], ()), "ggt": ("BIOPRO", ["LBXSGTSI"], ()), "albumin_serum": ("BIOPRO", ["LBXSAL"], ()),
    "bun": ("BIOPRO", ["LBXSBU"], ()), "calcium": ("BIOPRO", ["LBXSCA"], ()),
    "u_alb": ("ALB_CR", ["URXUMA", "URXUMS"], ()), "u_cr": ("ALB_CR", ["URXUCR"], ()), "uacr_cdc": ("ALB_CR", ["URDACT"], ()),
    "wbc": ("CBC", ["LBXWBCSI"], ()), "lymph_pct": ("CBC", ["LBXLYPCT"], ()), "mono_pct": ("CBC", ["LBXMOPCT"], ()),
    "neut_pct": ("CBC", ["LBXNEPCT"], ()), "eos_pct": ("CBC", ["LBXEOPCT"], ()), "baso_pct": ("CBC", ["LBXBAPCT"], ()),
    "lymph_n": ("CBC", ["LBDLYMNO"], ()), "mono_n": ("CBC", ["LBDMONO"], ()), "neut_n": ("CBC", ["LBDNENO"], ()),
    "eos_n": ("CBC", ["LBDEONO"], ()), "baso_n": ("CBC", ["LBDBANO"], ()), "rbc": ("CBC", ["LBXRBCSI"], ()),
    "hgb": ("CBC", ["LBXHGB"], ()), "hct": ("CBC", ["LBXHCT"], ()), "mcv": ("CBC", ["LBXMCVSI"], ()),
    "rdw": ("CBC", ["LBXRDW"], ()), "platelets": ("CBC", ["LBXPLTSI"], ()),
}
BP_AUSC = ("BPX", ["BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"], ["BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4"])
BP_OSC = ("BPXO", ["BPXOSY1", "BPXOSY2", "BPXOSY3"], ["BPXODI1", "BPXODI2", "BPXODI3"])

def source_of(cycle, name):
    """Return 'FILE.COLUMN' for a harmonised variable in a cycle, or None if absent."""
    comp, cols, _ = VARS[name]
    x = read(cycle, comp)
    if x is None: return None
    for c in cols:
        if c in x.columns: return f"{os.path.basename(file_path(cycle, comp))}:{c}"
    return None

# ---------------------------------------------------------------- equations
SCR_CAL = {"1999-2000": (1.013, 0.147), "2005-2006": (0.978, -0.016)}  # (slope, intercept)

def egfr_ckdepi2009(scr, age, female, black):
    female = np.asarray(female, bool); k = np.where(female, 0.7, 0.9); a = np.where(female, -0.329, -0.411); r = np.asarray(scr, float) / k
    return 141 * np.minimum(r, 1) ** a * np.maximum(r, 1) ** -1.209 * 0.993 ** np.asarray(age, float) * np.where(female, 1.018, 1.0) * np.where(np.asarray(black, bool), 1.159, 1.0)

def egfr_ckdepi2021(scr, age, female):
    female = np.asarray(female, bool); k = np.where(female, 0.7, 0.9); a = np.where(female, -0.241, -0.302); r = np.asarray(scr, float) / k
    return 142 * np.minimum(r, 1) ** a * np.maximum(r, 1) ** -1.200 * 0.9938 ** np.asarray(age, float) * np.where(female, 1.012, 1.0)

def _flag(cond, observed):
    return pd.Series(np.where(cond, 1.0, np.where(observed, 0.0, np.nan)), index=cond.index)

# ---------------------------------------------------------------- builder
def load_cycle(cycle):
    demo = read(cycle, "DEMO")
    if demo is None: raise FileNotFoundError(f"DEMO missing for {cycle}")
    df = pd.DataFrame({"SEQN": demo.SEQN.values}); df["cycle"] = cycle
    src = {}
    for name, (comp, cols, codes) in VARS.items():
        x = read(cycle, comp); col = None
        if x is not None: col = next((c for c in cols if c in x.columns), None)
        if col is None: df[name] = np.nan; continue
        s = x.drop_duplicates("SEQN").set_index("SEQN")[col].reindex(df.SEQN).to_numpy(float)
        if codes: s = np.where(np.isin(s, codes), np.nan, s)
        df[name] = s; src[name] = col
    # blood pressure
    osc = CYCLES[cycle][5]
    for tag, (comp, sy, di) in [("ausc", BP_AUSC), ("osc", BP_OSC)]:
        x = read(cycle, comp)
        if x is None: df[f"sbp_{tag}"] = np.nan; df[f"dbp_{tag}"] = np.nan; continue
        x = x.drop_duplicates("SEQN").set_index("SEQN").reindex(df.SEQN)
        sy_ = [c for c in sy if c in x.columns]; di_ = [c for c in di if c in x.columns]
        df[f"sbp_{tag}"] = x[sy_].mean(axis=1).to_numpy() if sy_ else np.nan
        df[f"dbp_{tag}"] = x[di_].replace(0, np.nan).mean(axis=1).to_numpy() if di_ else np.nan
    df["sbp"] = df.sbp_osc if osc else df.sbp_ausc
    df["dbp"] = df.dbp_osc if osc else df.dbp_ausc
    df["bp_protocol"] = "oscillometric" if osc else "auscultatory"
    # insulin: separate INS file from 2013, inside the fasting glucose file before
    df["insulin"] = df.insulin.fillna(df.pop("insulin_glu"))
    df["fast_hours"] = df.fast_hours.fillna(df.pop("fast_hours_fq"))   # FASTQX file (downloaded for 2021-2023 only)
    # creatinine calibration
    b, a = SCR_CAL.get(cycle, (1.0, 0.0))
    df["scr"] = a + b * df.scr_raw
    df["scr_calibrated"] = cycle in SCR_CAL
    return df

def derive(df):
    df = df.copy()
    df["female"] = (df.sex == 2).astype(float).where(df.sex.notna())
    df["pregnant"] = (df.pregnant_raw == 1).astype(float)
    df["race1"] = df.RIDRETH1.map({1: "mexican_american", 2: "other_hispanic", 3: "nh_white", 4: "nh_black", 5: "other_multi"})
    df["race3"] = df.RIDRETH3.map({1: "mexican_american", 2: "other_hispanic", 3: "nh_white", 4: "nh_black", 6: "nh_asian", 7: "other_multi"})
    df["educ3"] = df.educ_raw.map({1: "lt_hs", 2: "lt_hs", 3: "hs_ged", 4: "gt_hs", 5: "gt_hs"})
    df["smoking"] = np.select([df.SMQ020 == 2, (df.SMQ020 == 1) & (df.SMQ040 == 3), (df.SMQ020 == 1) & df.SMQ040.isin([1, 2])],
                              ["never", "former", "current"], default=None)
    df["smoking"] = df.smoking.where(df.smoking.notna(), np.nan)
    # diagnosed conditions (self-report); borderline diabetes (DIQ010 = 3) counts as not diagnosed
    df["dx_diabetes"] = _flag(df.DIQ010 == 1, df.DIQ010.isin([1, 2, 3]))
    df["dx_hypertension"] = _flag(df.BPQ020 == 1, df.BPQ020.isin([1, 2]))
    df["cvd_history"] = _flag((df[["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]] == 1).any(axis=1),
                              df[["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]].isin([1, 2]).all(axis=1))
    # composite flags kept for convenience (definitions vary by paper; papers' own definitions are rebuilt per paper)
    pos = (df.DIQ010 == 1) | (df.DIQ050 == 1) | (df.DIQ070 == 1) | (df.hba1c >= 6.5) | (df.glucose >= 126)
    df["diabetes_any"] = _flag(pos, df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna())
    for thr, (s, d) in {"140": (140, 90), "130": (130, 80)}.items():
        pos = (df.BPQ020 == 1) | (df.BPQ_med_now == 1) | (df.sbp >= s) | (df.dbp >= d)
        df[f"hypertension_any_{thr}"] = _flag(pos, df.BPQ020.isin([1, 2]) | (df.sbp.notna() & df.dbp.notna()))
    # kidney
    black = df.RIDRETH1 == 4
    df["egfr_ckdepi2009"] = egfr_ckdepi2009(df.scr, df.age, df.female == 1, black)
    df["egfr_ckdepi2021"] = egfr_ckdepi2021(df.scr, df.age, df.female == 1)
    df["uacr"] = np.where(df.u_cr > 0, df.u_alb / df.u_cr * 100, np.nan)
    ok = df.uacr.notna()
    df["albuminuria"] = (df.uacr >= 30).astype(float).where(ok)
    for eq in ["2009", "2021"]:
        e = df[f"egfr_ckdepi{eq}"]
        df[f"egfr_lt60_{eq}"] = (e < 60).astype(float).where(e.notna())
        df[f"ckd_{eq}"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(ok & e.notna())
    return df

def build(cycles=None):
    cycles = cycles or ALL_CYCLES
    if "2017-2018" in cycles and "2017-2020" in cycles:
        raise ValueError("2017-2018 (J) is contained in 2017-March 2020 (P); do not pool them")
    return derive(pd.concat([load_cycle(c) for c in cycles], ignore_index=True))

# ---------------------------------------------------------------- weights
WKIND = {"interview": ("WTINT2YR", "WTINT4YR", "WTINTPRP"), "mec": ("WTMEC2YR", "WTMEC4YR", "WTMECPRP"),
         "fasting": ("WTSAF2YR", "WTSAF4YR", "WTSAFPRP"),
         # 2021-2023 only: blood-draw components (BIOPRO, GHB, HDL, TCHOL, CBC) carry the phlebotomy weight WTPH2YR
         "phlebotomy": ("WTPH2YR", None, None)}

def pooled_weight(df, kind):
    """Pooled weight for the cycles present in df (NCHS rule; see module docstring). Rows outside the subsample get 0."""
    w2, w4, wp = WKIND[kind]
    cyc = list(dict.fromkeys(df.cycle))
    if "2017-2018" in cyc and "2017-2020" in cyc: raise ValueError("J and P overlap")
    four = "1999-2000" in cyc and "2001-2002" in cyc
    total = sum(CYCLES[c][4] for c in cyc)
    w = pd.Series(0.0, index=df.index)
    for c in cyc:
        m = df.cycle == c
        if c == "2017-2020": col, f = wp, 3.2 / total
        elif four and c in ("1999-2000", "2001-2002"): col, f = w4, 4 / total   # WT*4YR spans 1999-2002 (4 years)
        else: col, f = w2, 2 / total
        if kind == "phlebotomy" and c != "2021-2023":   # before 2021 blood-draw components use the MEC weight
            col = "WTMECPRP" if c == "2017-2020" else ("WTMEC4YR" if four and c in ("1999-2000", "2001-2002") else "WTMEC2YR")
        w[m] = df.loc[m, col].fillna(0) * f
    return w

def check_strata(df, w):
    """Strata in the pooled design with fewer than 2 PSUs among rows with positive weight (must be empty)."""
    d = df.loc[w > 0, ["SDMVSTRA", "SDMVPSU"]].drop_duplicates()
    n = d.groupby("SDMVSTRA").size()
    return n[n < 2]

# ---------------------------------------------------------------- estimation
def _irls(X, y, w, it=50):
    beta = np.zeros(X.shape[1])
    for _ in range(it):
        p = 1 / (1 + np.exp(-np.clip(X @ beta, -30, 30)))
        H = X.T @ (X * (w * p * (1 - p))[:, None]); g = X.T @ (w * (y - p))
        step = np.linalg.solve(H, g); beta += step
        if np.max(np.abs(step)) < 1e-10: break
    p = 1 / (1 + np.exp(-np.clip(X @ beta, -30, 30)))
    H = X.T @ (X * (w * p * (1 - p))[:, None])
    return beta, p, H

def svy_logit(df, y, x, w, domain=None, strata="SDMVSTRA", psu="SDMVPSU", df_rule="design"):
    """Weighted logistic regression with Taylor-linearisation (stratified, with-replacement PSU) variance.
    df: ALL rows of the pooled cycles (needed so that every PSU of the design is kept); domain: bool mask of the
    subpopulation; complete cases on y and x within the domain are analysed (domain estimation: other rows get zero score).
    df_rule: "design" = #PSU - #strata (prereg section 6), "svyglm" = design df + 1 - #parameters (sensitivity).
    Returns dict with coef, se, names, n, df, and a summary DataFrame (OR, t-based 95% CI)."""
    wt = (df[w] if isinstance(w, str) else pd.Series(w, index=df.index)).fillna(0).to_numpy(float)
    dom = np.ones(len(df), bool) if domain is None else np.asarray(domain, bool)
    cc = dom & df[[y] + list(x)].notna().all(axis=1).to_numpy() & (wt > 0)
    S = df[cc]; X = np.column_stack([np.ones(cc.sum()), S[list(x)].to_numpy(float)]); Y = S[y].to_numpy(float); W = wt[cc] / wt[cc].mean()
    beta, p, H = _irls(X, Y, W)
    U = X * (W * (Y - p))[:, None]
    des = df[[strata, psu]].drop_duplicates()   # every PSU of the pooled design (domain estimation)
    key = lambda a, b: a.astype(np.int64) * 1000 + b.astype(np.int64)
    Ud = pd.DataFrame(U).groupby(key(S[strata].to_numpy(), S[psu].to_numpy())).sum()
    allk = key(des[strata].to_numpy(), des[psu].to_numpy())
    Ud = Ud.reindex(allk).fillna(0.0).to_numpy(); hs = des[strata].to_numpy()
    M = np.zeros((X.shape[1], X.shape[1]))
    for h in np.unique(hs):
        G = Ud[hs == h]; n = len(G)
        if n < 2: continue
        Gc = G - G.mean(0); M += n / (n - 1) * Gc.T @ Gc
    Hi = np.linalg.inv(H); V = Hi @ M @ Hi; se = np.sqrt(np.diag(V))
    ddf = len(allk) - len(np.unique(hs))
    if df_rule == "svyglm": ddf = ddf + 1 - X.shape[1]
    q = stats.t.ppf(0.975, ddf); names = ["(Intercept)"] + list(x)
    tab = pd.DataFrame({"coef": beta, "se": se, "OR": np.exp(beta), "lo": np.exp(beta - q * se), "hi": np.exp(beta + q * se),
                        "p": 2 * stats.t.sf(np.abs(beta / se), ddf)}, index=names)
    return {"coef": beta, "se": se, "V": V, "names": names, "n": int(cc.sum()), "df": ddf, "table": tab}

def svy_mean(df, y, w, domain=None, strata="SDMVSTRA", psu="SDMVPSU"):
    """Weighted mean / proportion of y in a domain with Taylor-linearisation SE (ratio estimator) and design df."""
    wt = (df[w] if isinstance(w, str) else pd.Series(w, index=df.index)).fillna(0).to_numpy(float)
    dom = np.ones(len(df), bool) if domain is None else np.asarray(domain, bool)
    cc = dom & df[y].notna().to_numpy() & (wt > 0)
    yy = np.where(cc, df[y].fillna(0).to_numpy(float), 0.0); ww = np.where(cc, wt, 0.0)
    est = (ww * yy).sum() / ww.sum()
    z = ww * (yy - est) / ww.sum()
    k = df[strata].to_numpy(np.int64) * 1000 + df[psu].to_numpy(np.int64)
    t = pd.DataFrame({"h": df[strata].to_numpy(), "k": k, "z": z}).groupby(["h", "k"]).z.sum().reset_index()
    var = 0.0
    for h, g in t.groupby("h"):
        n = len(g)
        if n > 1: var += n / (n - 1) * ((g.z - g.z.mean()) ** 2).sum()
    ddf = t.k.nunique() - t.h.nunique()
    return {"est": est, "se": np.sqrt(var), "n": int(cc.sum()), "df": ddf}

# ---------------------------------------------------------------- documentation
def write_doc(path=os.path.join(HERE, "nhanes_core_variables.md")):
    head = open(path, encoding="utf-8").read().split("<!-- AUTO -->")[0] if os.path.exists(path) else ""
    lines = ["<!-- AUTO -->", "", "## Source column per cycle (generated by `python -I nhanes_core.py --doc` from the XPT column catalogue)", "",
             "Cell = file:column used; `absent` = not in the downloaded CDC files for that cycle.", ""]
    cyc = ALL_CYCLES
    lines.append("| Variable | " + " | ".join(cyc) + " |"); lines.append("|---" * (len(cyc) + 1) + "|")
    for name in VARS:
        lines.append(f"| {name} | " + " | ".join(source_of(c, name) or "absent" for c in cyc) + " |")
    for tag, (comp, sy, di) in [("sbp/dbp auscultatory", BP_AUSC), ("sbp/dbp oscillometric", BP_OSC)]:
        cells = []
        for c in cyc:
            x = read(c, comp)
            cells.append("absent" if x is None else f"{os.path.basename(file_path(c, comp))}:" + ",".join(v for v in sy + di if v in x.columns))
        lines.append(f"| {tag} | " + " | ".join(cells) + " |")
    open(path, "w", encoding="utf-8").write(head + "\n".join(lines) + "\n")

if __name__ == "__main__":
    if "--doc" in sys.argv: write_doc(); print("doc written")
    if "--summary" in sys.argv:
        for c in ALL_CYCLES:
            d = build([c])
            print(c, len(d), "adults", int((d.age >= 20).sum()), "uacr", int(d.uacr.notna().sum()), "scr", int(d.scr.notna().sum()),
                  "fastw>0", int((pooled_weight(d, "fasting") > 0).sum()), "tg", int(d.tg.notna().sum()), "hba1c", int(d.hba1c.notna().sum()),
                  "sbp", int(d.sbp.notna().sum()), "wbc", int(d.wbc.notna().sum()), flush=True)

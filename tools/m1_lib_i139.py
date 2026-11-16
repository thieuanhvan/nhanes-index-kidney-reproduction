"""I139 M1 helper library: shared loading for the per-paper specification builders in tools/m1_specs/<rec_id>.py.

A builder module defines
    REC_ID = "<rec_id>"
    def build(lib):  -> dict with keys
        df        : DataFrame of ALL rows of the paper's pooled cycles (needed for domain estimation; do not subset rows)
        domain    : bool Series, the paper's analytic population before complete-case on y and model terms
        y         : name of the 0/1 outcome column
        terms     : list of model term columns (exposure terms first, then covariate columns; dummies already made)
        target    : name of the exposure term whose OR is compared with the paper (e.g. "x" or "x_q4")
        weight    : name of the weight column, or None when the paper reports no weighting (then w = 1)
        readings  : dict point -> reading used (reading_1 of the locked spec unless a new point arises; log it)
        notes     : short text
The builder must not fit any model and must not look at model output; scripts/m1_run_i139.py fits and compares.

lib functions:
    table(cycles)               harmonised nhanes_core table + height/weight + extra components (cached)
    raw(df, stem_by_cycle, col)  pull any column of any downloaded XPT, matched on cycle and SEQN
                                 stem_by_cycle: dict cycle -> file stem (e.g. {"2011-2012": "ALQ_G"}) or a suffix rule
                                 via lib.stems("ALQ") which gives the usual names for all cycles
    rx_any(df, pattern)         1/0: any prescription in RXQ_RX whose drug name (RXDDRUG) matches the regex pattern
                                 (NaN when the person has no RXQ_RX record)
    dummies(df, col, ref)       0/1 columns for each level except ref, NaN when col is NaN; returns list of names
    quantile_dummies(df, col, mask, k, weight=None) quantile group dummies (Q2..Qk vs Q1) cut on the domain rows
    sd_scale(df, col, mask)     column / SD over the domain rows (unweighted)
    INDEX, DIET                 index functions (recon_i139.INDEX keys, diet_indices_i139.DIET keys)
    pooled_weight(df, kind)     nhanes_core pooled weight ("interview", "mec", "fasting", "phlebotomy")
"""
import os, sys, re
import numpy as np, pandas as pd
# Fix 15/11/2026 (extraction/m1_correction_20261115.md): with numexpr on (pandas default, numexpr 2.14.1, 8 threads),
# DataFrame.mask on frames above 1e6 cells intermittently returns all zeros. The XPT zero fix in nhanes_core.read,
# nhanes_extra_i139.read_stem and _read below then zeroes a whole large file (DEMO, DR1TOT, ...) and one cycle drops
# out of the sample without an error. Plain numpy evaluation gives the same numbers and is stable.
pd.set_option("compute.use_numexpr", False)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nhanes_core as nc
import nhanes_indices as ni
import nhanes_extra_i139 as nx
from recon_i139 import INDEX
try:
    from diet_indices_i139 import DIET
except Exception:  # dietary module optional
    DIET = {}

SUF = nx.SUF
pooled_weight = nc.pooled_weight
svy_logit = nc.svy_logit
_tab = {}


def table(cycles):
    k = tuple(cycles)
    if k not in _tab:
        _tab[k] = nx.add_extra(ni.add_anthropometry(nc.build(list(cycles))))
    return _tab[k].copy()


def stems(base):
    """Usual CDC file stems of a questionnaire/lab component for every release."""
    d = {c: base + s for c, s in SUF.items()}
    d["2017-2020"] = "P_" + base
    return d


_raw = {}


def _read(cycle, stem):
    p = os.path.join(nc.ROOT, cycle, stem + ".xpt")
    if not os.path.exists(p):
        return None
    if p not in _raw:
        x = pd.read_sas(p, format="xport")
        x.columns = [c.upper() for c in x.columns]
        num = x.select_dtypes("number").columns
        x[num] = x[num].mask(x[num].abs() < 1e-70, 0.0)
        x["SEQN"] = x["SEQN"].astype(np.int64)
        _raw[p] = x
    return _raw[p]


def raw(df, stem_by_cycle, col, first=True):
    """Column `col` (or first present of a list) from the file of each cycle; one row per SEQN (first record)."""
    out = pd.Series(np.nan, index=df.index, dtype=float)
    cols = [col] if isinstance(col, str) else list(col)
    for c in df.cycle.unique():
        st = stem_by_cycle.get(c)
        x = None if st is None else _read(c, st)
        if x is None:
            continue
        cc = next((k for k in cols if k in x.columns), None)
        if cc is None:
            continue
        x = x.drop_duplicates("SEQN").set_index("SEQN")
        m = (df.cycle == c).to_numpy()
        out[m] = pd.to_numeric(x[cc], errors="coerce").reindex(df.loc[m, "SEQN"]).to_numpy(float)
    return out


def rx_any(df, pattern):
    out = pd.Series(np.nan, index=df.index, dtype=float)
    rx = re.compile(pattern, re.I)
    for c in df.cycle.unique():
        st = ("RXQ_RX" + SUF[c]) if c in SUF else "P_RXQ_RX"
        x = _read(c, st)
        if x is None or "RXDDRUG" not in x.columns:
            continue
        names = x["RXDDRUG"].astype(str)
        if names.str.startswith("b'").any():
            names = x["RXDDRUG"].map(lambda v: v.decode("latin-1") if isinstance(v, bytes) else str(v))
        hit = x.assign(h=names.str.contains(rx)).groupby("SEQN").h.any().astype(float)
        m = (df.cycle == c).to_numpy()
        out[m] = hit.reindex(df.loc[m, "SEQN"]).to_numpy(float)
    return out


def dummies(df, col, ref):
    names = []
    for lv in sorted(df[col].dropna().unique(), key=str):
        if lv == ref:
            continue
        n = f"{col}__{lv}"
        df[n] = (df[col] == lv).astype(float).where(df[col].notna())
        names.append(n)
    return names


def quantile_dummies(df, col, mask, k=4, weight=None, prefix=None):
    v = df.loc[mask & df[col].notna(), col]
    if weight is None:
        cuts = np.quantile(v, np.linspace(0, 1, k + 1)[1:-1])
    else:
        w = df.loc[v.index, weight].to_numpy(float)
        o = np.argsort(v.to_numpy())
        cw = np.cumsum(w[o]) / w.sum()
        cuts = [v.to_numpy()[o][np.searchsorted(cw, q)] for q in np.linspace(0, 1, k + 1)[1:-1]]
    g = pd.Series(np.digitize(df[col], cuts, right=True) + 1, index=df.index).where(df[col].notna())
    prefix = prefix or f"{col}_q"
    names = []
    for j in range(2, k + 1):
        df[f"{prefix}{j}"] = (g == j).astype(float).where(g.notna())
        names.append(f"{prefix}{j}")
    df[prefix + "group"] = g
    return names, list(cuts)


def sd_scale(df, col, mask):
    sd = df.loc[mask & df[col].notna(), col].std()
    return df[col] / sd

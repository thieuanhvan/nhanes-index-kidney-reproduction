"""Tests for tools/nhanes_indices.py: three synthetic persons per index, expected values computed by hand
(independent scalar arithmetic, not the module). Run: python -I -m pytest tests/test_nhanes_indices.py
Person 2 is built so that TG and HDL are both 1 mmol/L (AIP = 0, LAP = WC - 58, CMI = WHtR)."""
import math, os, sys
import numpy as np, pandas as pd, pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import nhanes_indices as ni

PERSONS = pd.DataFrame([
    dict(female=0, tg=150, glucose=100, hdl=40, tc=200, bmi=27, waist=95, height=175, weight=82.7, insulin=10, hba1c=5.6,
         hypertension_any_140=0, neut_n=4.0, lymph_n=2.0, mono_n=0.5, platelets=250, neut_pct=60, albumin_serum=4.2, hgb=15.0),
    dict(female=1, tg=88.57, glucose=90, hdl=38.67, tc=154.68, bmi=22, waist=70, height=160, weight=56.3, insulin=5, hba1c=5.0,
         hypertension_any_140=1, neut_n=3.0, lymph_n=3.0, mono_n=0.3, platelets=300, neut_pct=50, albumin_serum=4.5, hgb=13.0),
    dict(female=1, tg=300, glucose=140, hdl=50, tc=260, bmi=35, waist=112, height=165, weight=95.3, insulin=25, hba1c=7.2,
         hypertension_any_140=1, neut_n=6.5, lymph_n=1.3, mono_n=0.8, platelets=180, neut_pct=72, albumin_serum=3.6, hgb=11.5),
])

EXPECTED = {
    "TyG": (8.9226583, 8.29045569, 9.952277717),
    "TyG-BMI": (240.9117741, 182.3900252, 348.3297201),
    "TyG-WC": (847.6525385, 580.3318983, 1114.655104),
    "TyG-WHtR": (4.843728791, 3.627074364, 6.75548548),
    "NLR": (2, 1, 5),
    "PLR": (125, 100, 138.4615385),
    "MLR": (0.25, 0.1, 0.6153846154),
    "SII": (500, 300, 900),
    "SIRI": (1, 0.3, 4),
    "AISI": (250, 90, 720),
    "NPAR": (14.28571429, 11.11111111, 20),
    "HALP": (50.4, 58.5, 29.9),
    "AIP": (0.2141187951, 0, 0.4182387778),
    "NHHR": (4, 3, 4.2),
    "VAI": (2.1873381, 1.680629794, 5.359423564),
    "LAP": (50.80727109, 12, 182.9061759),
    "CMI": (0.8888006258, 0.4375, 1.77816827),
    "BRI": (4.198035474, 2.261515664, 7.333473513),
    "ABSI": (0.07979249986, 0.07048384619, 0.0814886934),
    "WWI": (10.44650101, 9.329187949, 11.47285679),
    "WHtR": (0.5428571429, 0.4375, 0.6787878788),
    "RFM": (27.15789474, 30.28571429, 46.53571429),
    "HOMA-IR": (2.469135802, 1.111111111, 8.641975309),
    "HOMA-B": (97.2972973, 66.66666667, 116.8831169),
    "QUICKI": (0.3333333333, 0.3769015843, 0.2821616254),
    "METS-IR": (42.8759457, 33.66519748, 56.92859764),
    "eGDR": (9.5224, 8.696, 3.7038),
}


def test_every_index_has_expected_values():
    assert set(EXPECTED) == set(ni.INDICES)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_index_three_persons(name):
    got = ni.INDICES[name][1](PERSONS).to_numpy(float)
    for g, e in zip(got, EXPECTED[name]):
        assert math.isclose(g, e, rel_tol=1e-8, abs_tol=1e-9), (name, g, e)


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_missing_input_gives_nan(name):
    d = PERSONS.copy().astype(float)
    d.loc[0, :] = np.nan
    out = ni.INDICES[name][1](d)
    assert np.isnan(out.iloc[0]) and out.iloc[1:].notna().all()


def test_outcome_component_flags():
    flagged = {k for k, (_, _, flag) in ni.INDICES.items() if flag}
    assert flagged == {"NPAR", "HALP"}
    for k, (_, f, flag) in ni.INDICES.items():
        assert f"OUTCOME_COMPONENT: {flag}" in f.__doc__, k
        assert "https://doi.org/10." in f.__doc__, k


def test_guards():
    d = PERSONS.copy().astype(float)
    d.loc[0, "glucose"] = 60          # HOMA-B undefined at glucose <= 63 mg/dL
    d.loc[1, "lymph_n"] = 0           # zero denominator
    assert np.isnan(ni.homa_b(d).iloc[0])
    assert np.isnan(ni.nlr(d).iloc[1])


@pytest.mark.parametrize("f", [ni.life_simple_7, ni.life_essential_8])
def test_d3_stubs(f):
    with pytest.raises(NotImplementedError):
        f(PERSONS)


def test_trunc4():
    assert ni.trunc4(8.68599) == "8.685" and ni.trunc4(0.0799999) == "0.07999" and ni.trunc4(1327.9) == "1327"
    assert ni.trunc4(-2.34567) == "-2.345" and ni.trunc4(float("nan")) == "NA"

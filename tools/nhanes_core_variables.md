# nhanes_core: harmonised NHANES variables for I139

Module: `tools/nhanes_core.py`. Data: `data-external/nhanes-extra/<cycle>/` (public CDC/NCHS XPT files; URL, date and sha256 of each file in `download_log.json`). Pre-registration: `prereg_i139_v7_en.md`, sections 4 and 6.

## Releases covered

| Release | File naming | Years for pooling | Blood pressure protocol |
|---|---|---|---|
| 1999-2000 | DEMO, LAB16, LAB18, LAB10AM, LAB13AM, LAB13, LAB10, LAB25, BMX, BPX, DIQ, BPQ, SMQ, MCQ | 2 (4 with 2001-2002, WT*4YR) | auscultatory |
| 2001-2002 | L* names with suffix _B (L16, L40, L10AM, L13AM, L13, L10, L25) plus DEMO_B, BMX_B, BPX_B, DIQ_B, BPQ_B, SMQ_B, MCQ_B | 2 (4 with 1999-2000) | auscultatory |
| 2003-2004 | as 2001-2002, suffix _C | 2 | auscultatory |
| 2005-2006 to 2015-2016 | ALB_CR, BIOPRO, GLU, TRIGLY, HDL, TCHOL, GHB, CBC (+ INS from 2013), suffix _D to _I | 2 | auscultatory |
| 2017-2018 | suffix _J; BPX_J and BPXO_J both present | 2 | oscillometric (BPXO_J); auscultatory kept as sbp_ausc |
| 2017-March 2020 pre-pandemic | prefix P_ | 3.2 (WT*PRP) | oscillometric |
| 2021-2023 | suffix _L | 2 | oscillometric |

2017-2018 (J) is a subset of 2017-March 2020 (P): `build()` and `pooled_weight()` refuse to pool them.

## Derived variables and rules

- `female`, `pregnant` (RIDEXPRG = 1), `race1` (RIDRETH1, 5 groups), `race3` (RIDRETH3, 6 groups, 2011 onwards), `educ3` (DMDEDUC2: 1-2 less than high school, 3 high school or GED, 4-5 more than high school), `pir` (INDFMPIR).
- `smoking`: never (SMQ020 = 2), former (SMQ020 = 1 and SMQ040 = 3), current (SMQ020 = 1 and SMQ040 = 1 or 2).
- Missing codes 7, 9, 77, 99, 777, 999 (7777, 9999) of questionnaire items and codes 7, 9 of DMDEDUC2 are set to missing. No "missing" category is created.
- `dx_diabetes` = DIQ010 = 1 (borderline, code 3, counts as not diagnosed); `dx_hypertension` = BPQ020 = 1; `cvd_history` from MCQ160B to MCQ160F.
- Convenience composites (each paper's own definition is rebuilt per paper): `diabetes_any` (diagnosed, insulin, pills, HbA1c 6.5 % or more, fasting glucose 126 mg/dL or more); `hypertension_any_140` and `hypertension_any_130` (diagnosed, current medication BPQ050A, or BPQ150 in 2021-2023, or mean BP at the threshold).
- Blood pressure: `sbp`, `dbp` = mean of available readings of the protocol of the release (auscultatory BPXSY1-4/BPXDI1-4 before 2017, oscillometric BPXOSY1-3/BPXODI1-3 from 2017); diastolic 0 is set to missing. Both protocols are kept as `sbp_ausc`/`sbp_osc` where present.
- Serum creatinine `scr`: calibrated to the standard as recommended by NHANES: 1999-2000 `1.013 * Scr + 0.147`; 2005-2006 `-0.016 + 0.978 * Scr`. Other releases are used as published (the 2001-2002 column is LBDSCR). `scr_raw` keeps the published value, `scr_calibrated` flags calibrated rows.
- eGFR: `egfr_ckdepi2009` (race term for non-Hispanic Black, RIDRETH1 = 4) and `egfr_ckdepi2021` (race-free).
- `uacr` (mg/g) = urine albumin (mg/L) / urine creatinine (mg/dL) * 100; `uacr_cdc` keeps URDACT where CDC publishes it.
- Kidney flags: `albuminuria` (UACR 30 or more), `egfr_lt60_2009`, `egfr_lt60_2021`, `ckd_2009`, `ckd_2021` (UACR 30 or more, or eGFR below 60; missing when either input is missing).
- Insulin: INS file from 2013-2014, inside the fasting glucose file before (row `insulin_glu`). No cross-assay insulin calibration is applied (NHANES documents assay changes across 1999-2012); a paper that uses insulin is reconstructed with the paper's own handling and the assay change is noted.
- Fasting glucose and HbA1c are used as published (no cross-cycle calibration).
- Fasting hours `fast_hours`: PHAFSTHR from the GLU file (2005-2014) or FASTQX_L (2021-2023). FASTQX files of the other releases were not downloaded.

## Survey design and weights

- Design: `SDMVSTRA`, `SDMVPSU` (pseudo-PSU key = stratum * 1000 + PSU). `check_strata()` lists strata with fewer than 2 PSUs among rows with positive weight (must be empty before the Rao-Wu bootstrap; empty for 1999-2004 and for 2015-March 2020 in the test runs).
- `pooled_weight(df, kind)`, kind = `interview` (WTINT*), `mec` (WTMEC*), `fasting` (WTSAF*), `phlebotomy` (WTPH2YR in 2021-2023, where BIOPRO, GHB, HDL, TCHOL and CBC carry the phlebotomy weight; MEC weight in earlier releases). Rule: weight * years / total years of the pooled releases, with years = 2, or 3.2 for 2017-March 2020 (WT*PRP); for 1999-2000 and 2001-2002 pooled together, WT*4YR * 4 / total. Check: pooled MEC weights sum to 281.1 million for 1999-2004 and to 320.0 million for 2015-March 2020.
- In 2021-2023 the urine albumin/creatinine file has no own weight (MEC weight applies) and the fasting weight is in GLU_L, TRIGLY_L and INS_L.
- XPT zeros: pandas reads a stored 0 as 5.4e-79. `read()` sets these to exact 0, so "weight > 0" counts are correct (effect on pilot N: `pilot/kd3_variance_check_en.md`).
- Estimators: `svy_logit()` (weighted logistic, Taylor linearisation with strata and with-replacement PSUs, domain estimation keeping every PSU of the pooled design, df = PSUs minus strata, option `df_rule="svyglm"`), `svy_mean()` (weighted proportion with linearised SE). Checked in `pilot/kd3_variance_check_en.md`.

<!-- AUTO -->

## Source column per cycle (generated by `python -I nhanes_core.py --doc` from the XPT column catalogue)

Cell = file:column used; `absent` = not in the downloaded CDC files for that cycle.

| Variable | 1999-2000 | 2001-2002 | 2003-2004 | 2005-2006 | 2007-2008 | 2009-2010 | 2011-2012 | 2013-2014 | 2015-2016 | 2017-2018 | 2017-2020 | 2021-2023 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| age | DEMO.xpt:RIDAGEYR | DEMO_B.xpt:RIDAGEYR | DEMO_C.xpt:RIDAGEYR | DEMO_D.xpt:RIDAGEYR | DEMO_E.xpt:RIDAGEYR | DEMO_F.xpt:RIDAGEYR | DEMO_G.xpt:RIDAGEYR | DEMO_H.xpt:RIDAGEYR | DEMO_I.xpt:RIDAGEYR | DEMO_J.xpt:RIDAGEYR | P_DEMO.xpt:RIDAGEYR | DEMO_L.xpt:RIDAGEYR |
| sex | DEMO.xpt:RIAGENDR | DEMO_B.xpt:RIAGENDR | DEMO_C.xpt:RIAGENDR | DEMO_D.xpt:RIAGENDR | DEMO_E.xpt:RIAGENDR | DEMO_F.xpt:RIAGENDR | DEMO_G.xpt:RIAGENDR | DEMO_H.xpt:RIAGENDR | DEMO_I.xpt:RIAGENDR | DEMO_J.xpt:RIAGENDR | P_DEMO.xpt:RIAGENDR | DEMO_L.xpt:RIAGENDR |
| RIDRETH1 | DEMO.xpt:RIDRETH1 | DEMO_B.xpt:RIDRETH1 | DEMO_C.xpt:RIDRETH1 | DEMO_D.xpt:RIDRETH1 | DEMO_E.xpt:RIDRETH1 | DEMO_F.xpt:RIDRETH1 | DEMO_G.xpt:RIDRETH1 | DEMO_H.xpt:RIDRETH1 | DEMO_I.xpt:RIDRETH1 | DEMO_J.xpt:RIDRETH1 | P_DEMO.xpt:RIDRETH1 | DEMO_L.xpt:RIDRETH1 |
| RIDRETH3 | absent | absent | absent | absent | absent | absent | DEMO_G.xpt:RIDRETH3 | DEMO_H.xpt:RIDRETH3 | DEMO_I.xpt:RIDRETH3 | DEMO_J.xpt:RIDRETH3 | P_DEMO.xpt:RIDRETH3 | DEMO_L.xpt:RIDRETH3 |
| educ_raw | DEMO.xpt:DMDEDUC2 | DEMO_B.xpt:DMDEDUC2 | DEMO_C.xpt:DMDEDUC2 | DEMO_D.xpt:DMDEDUC2 | DEMO_E.xpt:DMDEDUC2 | DEMO_F.xpt:DMDEDUC2 | DEMO_G.xpt:DMDEDUC2 | DEMO_H.xpt:DMDEDUC2 | DEMO_I.xpt:DMDEDUC2 | DEMO_J.xpt:DMDEDUC2 | P_DEMO.xpt:DMDEDUC2 | DEMO_L.xpt:DMDEDUC2 |
| pir | DEMO.xpt:INDFMPIR | DEMO_B.xpt:INDFMPIR | DEMO_C.xpt:INDFMPIR | DEMO_D.xpt:INDFMPIR | DEMO_E.xpt:INDFMPIR | DEMO_F.xpt:INDFMPIR | DEMO_G.xpt:INDFMPIR | DEMO_H.xpt:INDFMPIR | DEMO_I.xpt:INDFMPIR | DEMO_J.xpt:INDFMPIR | P_DEMO.xpt:INDFMPIR | DEMO_L.xpt:INDFMPIR |
| pregnant_raw | DEMO.xpt:RIDEXPRG | DEMO_B.xpt:RIDEXPRG | DEMO_C.xpt:RIDEXPRG | DEMO_D.xpt:RIDEXPRG | DEMO_E.xpt:RIDEXPRG | DEMO_F.xpt:RIDEXPRG | DEMO_G.xpt:RIDEXPRG | DEMO_H.xpt:RIDEXPRG | DEMO_I.xpt:RIDEXPRG | DEMO_J.xpt:RIDEXPRG | P_DEMO.xpt:RIDEXPRG | DEMO_L.xpt:RIDEXPRG |
| SDMVSTRA | DEMO.xpt:SDMVSTRA | DEMO_B.xpt:SDMVSTRA | DEMO_C.xpt:SDMVSTRA | DEMO_D.xpt:SDMVSTRA | DEMO_E.xpt:SDMVSTRA | DEMO_F.xpt:SDMVSTRA | DEMO_G.xpt:SDMVSTRA | DEMO_H.xpt:SDMVSTRA | DEMO_I.xpt:SDMVSTRA | DEMO_J.xpt:SDMVSTRA | P_DEMO.xpt:SDMVSTRA | DEMO_L.xpt:SDMVSTRA |
| SDMVPSU | DEMO.xpt:SDMVPSU | DEMO_B.xpt:SDMVPSU | DEMO_C.xpt:SDMVPSU | DEMO_D.xpt:SDMVPSU | DEMO_E.xpt:SDMVPSU | DEMO_F.xpt:SDMVPSU | DEMO_G.xpt:SDMVPSU | DEMO_H.xpt:SDMVPSU | DEMO_I.xpt:SDMVPSU | DEMO_J.xpt:SDMVPSU | P_DEMO.xpt:SDMVPSU | DEMO_L.xpt:SDMVPSU |
| WTINT2YR | DEMO.xpt:WTINT2YR | DEMO_B.xpt:WTINT2YR | DEMO_C.xpt:WTINT2YR | DEMO_D.xpt:WTINT2YR | DEMO_E.xpt:WTINT2YR | DEMO_F.xpt:WTINT2YR | DEMO_G.xpt:WTINT2YR | DEMO_H.xpt:WTINT2YR | DEMO_I.xpt:WTINT2YR | DEMO_J.xpt:WTINT2YR | absent | DEMO_L.xpt:WTINT2YR |
| WTMEC2YR | DEMO.xpt:WTMEC2YR | DEMO_B.xpt:WTMEC2YR | DEMO_C.xpt:WTMEC2YR | DEMO_D.xpt:WTMEC2YR | DEMO_E.xpt:WTMEC2YR | DEMO_F.xpt:WTMEC2YR | DEMO_G.xpt:WTMEC2YR | DEMO_H.xpt:WTMEC2YR | DEMO_I.xpt:WTMEC2YR | DEMO_J.xpt:WTMEC2YR | absent | DEMO_L.xpt:WTMEC2YR |
| WTINT4YR | DEMO.xpt:WTINT4YR | DEMO_B.xpt:WTINT4YR | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent |
| WTMEC4YR | DEMO.xpt:WTMEC4YR | DEMO_B.xpt:WTMEC4YR | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent |
| WTINTPRP | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent | P_DEMO.xpt:WTINTPRP | absent |
| WTMECPRP | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent | P_DEMO.xpt:WTMECPRP | absent |
| SMQ020 | SMQ.xpt:SMQ020 | SMQ_B.xpt:SMQ020 | SMQ_C.xpt:SMQ020 | SMQ_D.xpt:SMQ020 | SMQ_E.xpt:SMQ020 | SMQ_F.xpt:SMQ020 | SMQ_G.xpt:SMQ020 | SMQ_H.xpt:SMQ020 | SMQ_I.xpt:SMQ020 | SMQ_J.xpt:SMQ020 | P_SMQ.xpt:SMQ020 | SMQ_L.xpt:SMQ020 |
| SMQ040 | SMQ.xpt:SMQ040 | SMQ_B.xpt:SMQ040 | SMQ_C.xpt:SMQ040 | SMQ_D.xpt:SMQ040 | SMQ_E.xpt:SMQ040 | SMQ_F.xpt:SMQ040 | SMQ_G.xpt:SMQ040 | SMQ_H.xpt:SMQ040 | SMQ_I.xpt:SMQ040 | SMQ_J.xpt:SMQ040 | P_SMQ.xpt:SMQ040 | SMQ_L.xpt:SMQ040 |
| DIQ010 | DIQ.xpt:DIQ010 | DIQ_B.xpt:DIQ010 | DIQ_C.xpt:DIQ010 | DIQ_D.xpt:DIQ010 | DIQ_E.xpt:DIQ010 | DIQ_F.xpt:DIQ010 | DIQ_G.xpt:DIQ010 | DIQ_H.xpt:DIQ010 | DIQ_I.xpt:DIQ010 | DIQ_J.xpt:DIQ010 | P_DIQ.xpt:DIQ010 | DIQ_L.xpt:DIQ010 |
| DIQ050 | DIQ.xpt:DIQ050 | DIQ_B.xpt:DIQ050 | DIQ_C.xpt:DIQ050 | DIQ_D.xpt:DIQ050 | DIQ_E.xpt:DIQ050 | DIQ_F.xpt:DIQ050 | DIQ_G.xpt:DIQ050 | DIQ_H.xpt:DIQ050 | DIQ_I.xpt:DIQ050 | DIQ_J.xpt:DIQ050 | P_DIQ.xpt:DIQ050 | DIQ_L.xpt:DIQ050 |
| DIQ070 | DIQ.xpt:DIQ070 | DIQ_B.xpt:DIQ070 | DIQ_C.xpt:DIQ070 | absent | absent | DIQ_F.xpt:DIQ070 | DIQ_G.xpt:DIQ070 | DIQ_H.xpt:DIQ070 | DIQ_I.xpt:DIQ070 | DIQ_J.xpt:DIQ070 | P_DIQ.xpt:DIQ070 | DIQ_L.xpt:DIQ070 |
| BPQ020 | BPQ.xpt:BPQ020 | BPQ_B.xpt:BPQ020 | BPQ_C.xpt:BPQ020 | BPQ_D.xpt:BPQ020 | BPQ_E.xpt:BPQ020 | BPQ_F.xpt:BPQ020 | BPQ_G.xpt:BPQ020 | BPQ_H.xpt:BPQ020 | BPQ_I.xpt:BPQ020 | BPQ_J.xpt:BPQ020 | P_BPQ.xpt:BPQ020 | BPQ_L.xpt:BPQ020 |
| BPQ_med_now | BPQ.xpt:BPQ050A | BPQ_B.xpt:BPQ050A | BPQ_C.xpt:BPQ050A | BPQ_D.xpt:BPQ050A | BPQ_E.xpt:BPQ050A | BPQ_F.xpt:BPQ050A | BPQ_G.xpt:BPQ050A | BPQ_H.xpt:BPQ050A | BPQ_I.xpt:BPQ050A | BPQ_J.xpt:BPQ050A | P_BPQ.xpt:BPQ050A | BPQ_L.xpt:BPQ150 |
| MCQ160B | MCQ.xpt:MCQ160B | MCQ_B.xpt:MCQ160B | MCQ_C.xpt:MCQ160B | MCQ_D.xpt:MCQ160B | MCQ_E.xpt:MCQ160B | MCQ_F.xpt:MCQ160B | MCQ_G.xpt:MCQ160B | MCQ_H.xpt:MCQ160B | MCQ_I.xpt:MCQ160B | MCQ_J.xpt:MCQ160B | P_MCQ.xpt:MCQ160B | MCQ_L.xpt:MCQ160B |
| MCQ160C | MCQ.xpt:MCQ160C | MCQ_B.xpt:MCQ160C | MCQ_C.xpt:MCQ160C | MCQ_D.xpt:MCQ160C | MCQ_E.xpt:MCQ160C | MCQ_F.xpt:MCQ160C | MCQ_G.xpt:MCQ160C | MCQ_H.xpt:MCQ160C | MCQ_I.xpt:MCQ160C | MCQ_J.xpt:MCQ160C | P_MCQ.xpt:MCQ160C | MCQ_L.xpt:MCQ160C |
| MCQ160D | MCQ.xpt:MCQ160D | MCQ_B.xpt:MCQ160D | MCQ_C.xpt:MCQ160D | MCQ_D.xpt:MCQ160D | MCQ_E.xpt:MCQ160D | MCQ_F.xpt:MCQ160D | MCQ_G.xpt:MCQ160D | MCQ_H.xpt:MCQ160D | MCQ_I.xpt:MCQ160D | MCQ_J.xpt:MCQ160D | P_MCQ.xpt:MCQ160D | MCQ_L.xpt:MCQ160D |
| MCQ160E | MCQ.xpt:MCQ160E | MCQ_B.xpt:MCQ160E | MCQ_C.xpt:MCQ160E | MCQ_D.xpt:MCQ160E | MCQ_E.xpt:MCQ160E | MCQ_F.xpt:MCQ160E | MCQ_G.xpt:MCQ160E | MCQ_H.xpt:MCQ160E | MCQ_I.xpt:MCQ160E | MCQ_J.xpt:MCQ160E | P_MCQ.xpt:MCQ160E | MCQ_L.xpt:MCQ160E |
| MCQ160F | MCQ.xpt:MCQ160F | MCQ_B.xpt:MCQ160F | MCQ_C.xpt:MCQ160F | MCQ_D.xpt:MCQ160F | MCQ_E.xpt:MCQ160F | MCQ_F.xpt:MCQ160F | MCQ_G.xpt:MCQ160F | MCQ_H.xpt:MCQ160F | MCQ_I.xpt:MCQ160F | MCQ_J.xpt:MCQ160F | P_MCQ.xpt:MCQ160F | MCQ_L.xpt:MCQ160F |
| bmi | BMX.xpt:BMXBMI | BMX_B.xpt:BMXBMI | BMX_C.xpt:BMXBMI | BMX_D.xpt:BMXBMI | BMX_E.xpt:BMXBMI | BMX_F.xpt:BMXBMI | BMX_G.xpt:BMXBMI | BMX_H.xpt:BMXBMI | BMX_I.xpt:BMXBMI | BMX_J.xpt:BMXBMI | P_BMX.xpt:BMXBMI | BMX_L.xpt:BMXBMI |
| waist | BMX.xpt:BMXWAIST | BMX_B.xpt:BMXWAIST | BMX_C.xpt:BMXWAIST | BMX_D.xpt:BMXWAIST | BMX_E.xpt:BMXWAIST | BMX_F.xpt:BMXWAIST | BMX_G.xpt:BMXWAIST | BMX_H.xpt:BMXWAIST | BMX_I.xpt:BMXWAIST | BMX_J.xpt:BMXWAIST | P_BMX.xpt:BMXWAIST | BMX_L.xpt:BMXWAIST |
| glucose | LAB10AM.xpt:LBXGLU | L10AM_B.xpt:LBXGLU | L10AM_C.xpt:LBXGLU | GLU_D.xpt:LBXGLU | GLU_E.xpt:LBXGLU | GLU_F.xpt:LBXGLU | GLU_G.xpt:LBXGLU | GLU_H.xpt:LBXGLU | GLU_I.xpt:LBXGLU | GLU_J.xpt:LBXGLU | P_GLU.xpt:LBXGLU | GLU_L.xpt:LBXGLU |
| WTSAF2YR | LAB10AM.xpt:WTSAF2YR | L10AM_B.xpt:WTSAF2YR | L10AM_C.xpt:WTSAF2YR | GLU_D.xpt:WTSAF2YR | GLU_E.xpt:WTSAF2YR | GLU_F.xpt:WTSAF2YR | GLU_G.xpt:WTSAF2YR | GLU_H.xpt:WTSAF2YR | GLU_I.xpt:WTSAF2YR | GLU_J.xpt:WTSAF2YR | absent | GLU_L.xpt:WTSAF2YR |
| WTSAF4YR | LAB10AM.xpt:WTSAF4YR | L10AM_B.xpt:WTSAF4YR | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent |
| WTSAFPRP | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent | P_GLU.xpt:WTSAFPRP | absent |
| fast_hours | absent | absent | absent | GLU_D.xpt:PHAFSTHR | GLU_E.xpt:PHAFSTHR | GLU_F.xpt:PHAFSTHR | GLU_G.xpt:PHAFSTHR | GLU_H.xpt:PHAFSTHR | absent | absent | absent | absent |
| fast_hours_fq | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent | FASTQX_L.xpt:PHAFSTHR |
| WTPH2YR | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent | absent | BIOPRO_L.xpt:WTPH2YR |
| insulin | absent | absent | absent | absent | absent | absent | absent | INS_H.xpt:LBXIN | INS_I.xpt:LBXIN | INS_J.xpt:LBXIN | P_INS.xpt:LBXIN | INS_L.xpt:LBXIN |
| insulin_glu | LAB10AM.xpt:LBXIN | L10AM_B.xpt:LBXIN | L10AM_C.xpt:LBXIN | GLU_D.xpt:LBXIN | GLU_E.xpt:LBXIN | GLU_F.xpt:LBXIN | GLU_G.xpt:LBXIN | absent | absent | absent | absent | absent |
| hba1c | LAB10.xpt:LBXGH | L10_B.xpt:LBXGH | L10_C.xpt:LBXGH | GHB_D.xpt:LBXGH | GHB_E.xpt:LBXGH | GHB_F.xpt:LBXGH | GHB_G.xpt:LBXGH | GHB_H.xpt:LBXGH | GHB_I.xpt:LBXGH | GHB_J.xpt:LBXGH | P_GHB.xpt:LBXGH | GHB_L.xpt:LBXGH |
| tg | LAB13AM.xpt:LBXTR | L13AM_B.xpt:LBXTR | L13AM_C.xpt:LBXTR | TRIGLY_D.xpt:LBXTR | TRIGLY_E.xpt:LBXTR | TRIGLY_F.xpt:LBXTR | TRIGLY_G.xpt:LBXTR | TRIGLY_H.xpt:LBXTR | TRIGLY_I.xpt:LBXTR | TRIGLY_J.xpt:LBXTR | P_TRIGLY.xpt:LBXTR | TRIGLY_L.xpt:LBXTLG |
| ldl | LAB13AM.xpt:LBDLDL | L13AM_B.xpt:LBDLDL | L13AM_C.xpt:LBDLDL | TRIGLY_D.xpt:LBDLDL | TRIGLY_E.xpt:LBDLDL | TRIGLY_F.xpt:LBDLDL | TRIGLY_G.xpt:LBDLDL | TRIGLY_H.xpt:LBDLDL | TRIGLY_I.xpt:LBDLDL | TRIGLY_J.xpt:LBDLDL | P_TRIGLY.xpt:LBDLDL | TRIGLY_L.xpt:LBDLDL |
| hdl | LAB13.xpt:LBDHDL | L13_B.xpt:LBDHDL | L13_C.xpt:LBXHDD | HDL_D.xpt:LBDHDD | HDL_E.xpt:LBDHDD | HDL_F.xpt:LBDHDD | HDL_G.xpt:LBDHDD | HDL_H.xpt:LBDHDD | HDL_I.xpt:LBDHDD | HDL_J.xpt:LBDHDD | P_HDL.xpt:LBDHDD | HDL_L.xpt:LBDHDD |
| tc | LAB13.xpt:LBXTC | L13_B.xpt:LBXTC | L13_C.xpt:LBXTC | TCHOL_D.xpt:LBXTC | TCHOL_E.xpt:LBXTC | TCHOL_F.xpt:LBXTC | TCHOL_G.xpt:LBXTC | TCHOL_H.xpt:LBXTC | TCHOL_I.xpt:LBXTC | TCHOL_J.xpt:LBXTC | P_TCHOL.xpt:LBXTC | TCHOL_L.xpt:LBXTC |
| scr_raw | LAB18.xpt:LBXSCR | L40_B.xpt:LBDSCR | L40_C.xpt:LBXSCR | BIOPRO_D.xpt:LBXSCR | BIOPRO_E.xpt:LBXSCR | BIOPRO_F.xpt:LBXSCR | BIOPRO_G.xpt:LBXSCR | BIOPRO_H.xpt:LBXSCR | BIOPRO_I.xpt:LBXSCR | BIOPRO_J.xpt:LBXSCR | P_BIOPRO.xpt:LBXSCR | BIOPRO_L.xpt:LBXSCR |
| tg_serum | LAB18.xpt:LBXSTR | L40_B.xpt:LBXSTR | L40_C.xpt:LBXSTR | BIOPRO_D.xpt:LBXSTR | BIOPRO_E.xpt:LBXSTR | BIOPRO_F.xpt:LBXSTR | BIOPRO_G.xpt:LBXSTR | BIOPRO_H.xpt:LBXSTR | BIOPRO_I.xpt:LBXSTR | BIOPRO_J.xpt:LBXSTR | P_BIOPRO.xpt:LBXSTR | BIOPRO_L.xpt:LBXSTR |
| glucose_serum | LAB18.xpt:LBXSGL | L40_B.xpt:LBXSGL | L40_C.xpt:LBXSGL | BIOPRO_D.xpt:LBXSGL | BIOPRO_E.xpt:LBXSGL | BIOPRO_F.xpt:LBXSGL | BIOPRO_G.xpt:LBXSGL | BIOPRO_H.xpt:LBXSGL | BIOPRO_I.xpt:LBXSGL | BIOPRO_J.xpt:LBXSGL | P_BIOPRO.xpt:LBXSGL | BIOPRO_L.xpt:LBXSGL |
| tc_serum | LAB18.xpt:LBXSCH | L40_B.xpt:LBXSCH | L40_C.xpt:LBXSCH | BIOPRO_D.xpt:LBXSCH | BIOPRO_E.xpt:LBXSCH | BIOPRO_F.xpt:LBXSCH | BIOPRO_G.xpt:LBXSCH | BIOPRO_H.xpt:LBXSCH | BIOPRO_I.xpt:LBXSCH | BIOPRO_J.xpt:LBXSCH | P_BIOPRO.xpt:LBXSCH | BIOPRO_L.xpt:LBXSCH |
| uric_acid | LAB18.xpt:LBXSUA | L40_B.xpt:LBXSUA | L40_C.xpt:LBXSUA | BIOPRO_D.xpt:LBXSUA | BIOPRO_E.xpt:LBXSUA | BIOPRO_F.xpt:LBXSUA | BIOPRO_G.xpt:LBXSUA | BIOPRO_H.xpt:LBXSUA | BIOPRO_I.xpt:LBXSUA | BIOPRO_J.xpt:LBXSUA | P_BIOPRO.xpt:LBXSUA | BIOPRO_L.xpt:LBXSUA |
| alt | LAB18.xpt:LBXSATSI | L40_B.xpt:LBXSATSI | L40_C.xpt:LBXSATSI | BIOPRO_D.xpt:LBXSATSI | BIOPRO_E.xpt:LBXSATSI | BIOPRO_F.xpt:LBXSATSI | BIOPRO_G.xpt:LBXSATSI | BIOPRO_H.xpt:LBXSATSI | BIOPRO_I.xpt:LBXSATSI | BIOPRO_J.xpt:LBXSATSI | P_BIOPRO.xpt:LBXSATSI | BIOPRO_L.xpt:LBXSATSI |
| ast | LAB18.xpt:LBXSASSI | L40_B.xpt:LBXSASSI | L40_C.xpt:LBXSASSI | BIOPRO_D.xpt:LBXSASSI | BIOPRO_E.xpt:LBXSASSI | BIOPRO_F.xpt:LBXSASSI | BIOPRO_G.xpt:LBXSASSI | BIOPRO_H.xpt:LBXSASSI | BIOPRO_I.xpt:LBXSASSI | BIOPRO_J.xpt:LBXSASSI | P_BIOPRO.xpt:LBXSASSI | BIOPRO_L.xpt:LBXSASSI |
| ggt | LAB18.xpt:LBXSGTSI | L40_B.xpt:LBXSGTSI | L40_C.xpt:LBXSGTSI | BIOPRO_D.xpt:LBXSGTSI | BIOPRO_E.xpt:LBXSGTSI | BIOPRO_F.xpt:LBXSGTSI | BIOPRO_G.xpt:LBXSGTSI | BIOPRO_H.xpt:LBXSGTSI | BIOPRO_I.xpt:LBXSGTSI | BIOPRO_J.xpt:LBXSGTSI | P_BIOPRO.xpt:LBXSGTSI | BIOPRO_L.xpt:LBXSGTSI |
| albumin_serum | LAB18.xpt:LBXSAL | L40_B.xpt:LBXSAL | L40_C.xpt:LBXSAL | BIOPRO_D.xpt:LBXSAL | BIOPRO_E.xpt:LBXSAL | BIOPRO_F.xpt:LBXSAL | BIOPRO_G.xpt:LBXSAL | BIOPRO_H.xpt:LBXSAL | BIOPRO_I.xpt:LBXSAL | BIOPRO_J.xpt:LBXSAL | P_BIOPRO.xpt:LBXSAL | BIOPRO_L.xpt:LBXSAL |
| bun | LAB18.xpt:LBXSBU | L40_B.xpt:LBXSBU | L40_C.xpt:LBXSBU | BIOPRO_D.xpt:LBXSBU | BIOPRO_E.xpt:LBXSBU | BIOPRO_F.xpt:LBXSBU | BIOPRO_G.xpt:LBXSBU | BIOPRO_H.xpt:LBXSBU | BIOPRO_I.xpt:LBXSBU | BIOPRO_J.xpt:LBXSBU | P_BIOPRO.xpt:LBXSBU | BIOPRO_L.xpt:LBXSBU |
| calcium | LAB18.xpt:LBXSCA | L40_B.xpt:LBXSCA | L40_C.xpt:LBXSCA | BIOPRO_D.xpt:LBXSCA | BIOPRO_E.xpt:LBXSCA | BIOPRO_F.xpt:LBXSCA | BIOPRO_G.xpt:LBXSCA | BIOPRO_H.xpt:LBXSCA | BIOPRO_I.xpt:LBXSCA | BIOPRO_J.xpt:LBXSCA | P_BIOPRO.xpt:LBXSCA | BIOPRO_L.xpt:LBXSCA |
| u_alb | LAB16.xpt:URXUMA | L16_B.xpt:URXUMA | L16_C.xpt:URXUMA | ALB_CR_D.xpt:URXUMA | ALB_CR_E.xpt:URXUMA | ALB_CR_F.xpt:URXUMA | ALB_CR_G.xpt:URXUMA | ALB_CR_H.xpt:URXUMA | ALB_CR_I.xpt:URXUMA | ALB_CR_J.xpt:URXUMA | P_ALB_CR.xpt:URXUMA | ALB_CR_L.xpt:URXUMA |
| u_cr | LAB16.xpt:URXUCR | L16_B.xpt:URXUCR | L16_C.xpt:URXUCR | ALB_CR_D.xpt:URXUCR | ALB_CR_E.xpt:URXUCR | ALB_CR_F.xpt:URXUCR | ALB_CR_G.xpt:URXUCR | ALB_CR_H.xpt:URXUCR | ALB_CR_I.xpt:URXUCR | ALB_CR_J.xpt:URXUCR | P_ALB_CR.xpt:URXUCR | ALB_CR_L.xpt:URXUCR |
| uacr_cdc | absent | absent | absent | absent | absent | ALB_CR_F.xpt:URDACT | ALB_CR_G.xpt:URDACT | ALB_CR_H.xpt:URDACT | ALB_CR_I.xpt:URDACT | ALB_CR_J.xpt:URDACT | P_ALB_CR.xpt:URDACT | ALB_CR_L.xpt:URDACT |
| wbc | LAB25.xpt:LBXWBCSI | L25_B.xpt:LBXWBCSI | L25_C.xpt:LBXWBCSI | CBC_D.xpt:LBXWBCSI | CBC_E.xpt:LBXWBCSI | CBC_F.xpt:LBXWBCSI | CBC_G.xpt:LBXWBCSI | CBC_H.xpt:LBXWBCSI | CBC_I.xpt:LBXWBCSI | CBC_J.xpt:LBXWBCSI | P_CBC.xpt:LBXWBCSI | CBC_L.xpt:LBXWBCSI |
| lymph_pct | LAB25.xpt:LBXLYPCT | L25_B.xpt:LBXLYPCT | L25_C.xpt:LBXLYPCT | CBC_D.xpt:LBXLYPCT | CBC_E.xpt:LBXLYPCT | CBC_F.xpt:LBXLYPCT | CBC_G.xpt:LBXLYPCT | CBC_H.xpt:LBXLYPCT | CBC_I.xpt:LBXLYPCT | CBC_J.xpt:LBXLYPCT | P_CBC.xpt:LBXLYPCT | CBC_L.xpt:LBXLYPCT |
| mono_pct | LAB25.xpt:LBXMOPCT | L25_B.xpt:LBXMOPCT | L25_C.xpt:LBXMOPCT | CBC_D.xpt:LBXMOPCT | CBC_E.xpt:LBXMOPCT | CBC_F.xpt:LBXMOPCT | CBC_G.xpt:LBXMOPCT | CBC_H.xpt:LBXMOPCT | CBC_I.xpt:LBXMOPCT | CBC_J.xpt:LBXMOPCT | P_CBC.xpt:LBXMOPCT | CBC_L.xpt:LBXMOPCT |
| neut_pct | LAB25.xpt:LBXNEPCT | L25_B.xpt:LBXNEPCT | L25_C.xpt:LBXNEPCT | CBC_D.xpt:LBXNEPCT | CBC_E.xpt:LBXNEPCT | CBC_F.xpt:LBXNEPCT | CBC_G.xpt:LBXNEPCT | CBC_H.xpt:LBXNEPCT | CBC_I.xpt:LBXNEPCT | CBC_J.xpt:LBXNEPCT | P_CBC.xpt:LBXNEPCT | CBC_L.xpt:LBXNEPCT |
| eos_pct | LAB25.xpt:LBXEOPCT | L25_B.xpt:LBXEOPCT | L25_C.xpt:LBXEOPCT | CBC_D.xpt:LBXEOPCT | CBC_E.xpt:LBXEOPCT | CBC_F.xpt:LBXEOPCT | CBC_G.xpt:LBXEOPCT | CBC_H.xpt:LBXEOPCT | CBC_I.xpt:LBXEOPCT | CBC_J.xpt:LBXEOPCT | P_CBC.xpt:LBXEOPCT | CBC_L.xpt:LBXEOPCT |
| baso_pct | LAB25.xpt:LBXBAPCT | L25_B.xpt:LBXBAPCT | L25_C.xpt:LBXBAPCT | CBC_D.xpt:LBXBAPCT | CBC_E.xpt:LBXBAPCT | CBC_F.xpt:LBXBAPCT | CBC_G.xpt:LBXBAPCT | CBC_H.xpt:LBXBAPCT | CBC_I.xpt:LBXBAPCT | CBC_J.xpt:LBXBAPCT | P_CBC.xpt:LBXBAPCT | CBC_L.xpt:LBXBAPCT |
| lymph_n | LAB25.xpt:LBDLYMNO | L25_B.xpt:LBDLYMNO | L25_C.xpt:LBDLYMNO | CBC_D.xpt:LBDLYMNO | CBC_E.xpt:LBDLYMNO | CBC_F.xpt:LBDLYMNO | CBC_G.xpt:LBDLYMNO | CBC_H.xpt:LBDLYMNO | CBC_I.xpt:LBDLYMNO | CBC_J.xpt:LBDLYMNO | P_CBC.xpt:LBDLYMNO | CBC_L.xpt:LBDLYMNO |
| mono_n | LAB25.xpt:LBDMONO | L25_B.xpt:LBDMONO | L25_C.xpt:LBDMONO | CBC_D.xpt:LBDMONO | CBC_E.xpt:LBDMONO | CBC_F.xpt:LBDMONO | CBC_G.xpt:LBDMONO | CBC_H.xpt:LBDMONO | CBC_I.xpt:LBDMONO | CBC_J.xpt:LBDMONO | P_CBC.xpt:LBDMONO | CBC_L.xpt:LBDMONO |
| neut_n | LAB25.xpt:LBDNENO | L25_B.xpt:LBDNENO | L25_C.xpt:LBDNENO | CBC_D.xpt:LBDNENO | CBC_E.xpt:LBDNENO | CBC_F.xpt:LBDNENO | CBC_G.xpt:LBDNENO | CBC_H.xpt:LBDNENO | CBC_I.xpt:LBDNENO | CBC_J.xpt:LBDNENO | P_CBC.xpt:LBDNENO | CBC_L.xpt:LBDNENO |
| eos_n | LAB25.xpt:LBDEONO | L25_B.xpt:LBDEONO | L25_C.xpt:LBDEONO | CBC_D.xpt:LBDEONO | CBC_E.xpt:LBDEONO | CBC_F.xpt:LBDEONO | CBC_G.xpt:LBDEONO | CBC_H.xpt:LBDEONO | CBC_I.xpt:LBDEONO | CBC_J.xpt:LBDEONO | P_CBC.xpt:LBDEONO | CBC_L.xpt:LBDEONO |
| baso_n | LAB25.xpt:LBDBANO | L25_B.xpt:LBDBANO | L25_C.xpt:LBDBANO | CBC_D.xpt:LBDBANO | CBC_E.xpt:LBDBANO | CBC_F.xpt:LBDBANO | CBC_G.xpt:LBDBANO | CBC_H.xpt:LBDBANO | CBC_I.xpt:LBDBANO | CBC_J.xpt:LBDBANO | P_CBC.xpt:LBDBANO | CBC_L.xpt:LBDBANO |
| rbc | LAB25.xpt:LBXRBCSI | L25_B.xpt:LBXRBCSI | L25_C.xpt:LBXRBCSI | CBC_D.xpt:LBXRBCSI | CBC_E.xpt:LBXRBCSI | CBC_F.xpt:LBXRBCSI | CBC_G.xpt:LBXRBCSI | CBC_H.xpt:LBXRBCSI | CBC_I.xpt:LBXRBCSI | CBC_J.xpt:LBXRBCSI | P_CBC.xpt:LBXRBCSI | CBC_L.xpt:LBXRBCSI |
| hgb | LAB25.xpt:LBXHGB | L25_B.xpt:LBXHGB | L25_C.xpt:LBXHGB | CBC_D.xpt:LBXHGB | CBC_E.xpt:LBXHGB | CBC_F.xpt:LBXHGB | CBC_G.xpt:LBXHGB | CBC_H.xpt:LBXHGB | CBC_I.xpt:LBXHGB | CBC_J.xpt:LBXHGB | P_CBC.xpt:LBXHGB | CBC_L.xpt:LBXHGB |
| hct | LAB25.xpt:LBXHCT | L25_B.xpt:LBXHCT | L25_C.xpt:LBXHCT | CBC_D.xpt:LBXHCT | CBC_E.xpt:LBXHCT | CBC_F.xpt:LBXHCT | CBC_G.xpt:LBXHCT | CBC_H.xpt:LBXHCT | CBC_I.xpt:LBXHCT | CBC_J.xpt:LBXHCT | P_CBC.xpt:LBXHCT | CBC_L.xpt:LBXHCT |
| mcv | LAB25.xpt:LBXMCVSI | L25_B.xpt:LBXMCVSI | L25_C.xpt:LBXMCVSI | CBC_D.xpt:LBXMCVSI | CBC_E.xpt:LBXMCVSI | CBC_F.xpt:LBXMCVSI | CBC_G.xpt:LBXMCVSI | CBC_H.xpt:LBXMCVSI | CBC_I.xpt:LBXMCVSI | CBC_J.xpt:LBXMCVSI | P_CBC.xpt:LBXMCVSI | CBC_L.xpt:LBXMCVSI |
| rdw | LAB25.xpt:LBXRDW | L25_B.xpt:LBXRDW | L25_C.xpt:LBXRDW | CBC_D.xpt:LBXRDW | CBC_E.xpt:LBXRDW | CBC_F.xpt:LBXRDW | CBC_G.xpt:LBXRDW | CBC_H.xpt:LBXRDW | CBC_I.xpt:LBXRDW | CBC_J.xpt:LBXRDW | P_CBC.xpt:LBXRDW | CBC_L.xpt:LBXRDW |
| platelets | LAB25.xpt:LBXPLTSI | L25_B.xpt:LBXPLTSI | L25_C.xpt:LBXPLTSI | CBC_D.xpt:LBXPLTSI | CBC_E.xpt:LBXPLTSI | CBC_F.xpt:LBXPLTSI | CBC_G.xpt:LBXPLTSI | CBC_H.xpt:LBXPLTSI | CBC_I.xpt:LBXPLTSI | CBC_J.xpt:LBXPLTSI | P_CBC.xpt:LBXPLTSI | CBC_L.xpt:LBXPLTSI |
| sbp/dbp auscultatory | BPX.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_B.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_C.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_D.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_E.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_F.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_G.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_H.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_I.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | BPX_J.xpt:BPXSY1,BPXSY2,BPXSY3,BPXSY4,BPXDI1,BPXDI2,BPXDI3,BPXDI4 | absent | absent |
| sbp/dbp oscillometric | absent | absent | absent | absent | absent | absent | absent | absent | absent | BPXO_J.xpt:BPXOSY1,BPXOSY2,BPXOSY3,BPXODI1,BPXODI2,BPXODI3 | P_BPXO.xpt:BPXOSY1,BPXOSY2,BPXOSY3,BPXODI1,BPXODI2,BPXODI3 | BPXO_L.xpt:BPXOSY1,BPXOSY2,BPXOSY3,BPXODI1,BPXODI2,BPXODI3 |

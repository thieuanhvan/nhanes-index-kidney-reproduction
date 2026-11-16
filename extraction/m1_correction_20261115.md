# M1 correction, 15 November 2026 (I139)

Two M1 rows were wrong: 9835a67923 and 69700c2b09. Each recorded run had silently lost one whole NHANES cycle (2007-2008). This file keeps the old rows, gives the cause, the evidence, the fix and the before/after numbers. Times are local (UTC+7).

## 1 Cause

`pandas` evaluates `DataFrame.mask` on frames of more than 1,000,000 cells with `numexpr` by default (here `numexpr` 2.14.1, 8 threads, `pandas` 2.3.3, Python 3.13.9). On this machine that evaluation sometimes returns a frame of all zeros. Every XPT loader in the pipeline applies `x[num] = x[num].mask(x[num].abs() < 1e-70, 0.0)` (the fix for XPT zeros read as 5.4e-79) to the whole file just after reading it: `nhanes_core.read`, `nhanes_extra_i139.read_stem` and `m1_lib_i139._read`. When the bad result hits one large file (DEMO, DR1TOT and DR2TOT are above the 1e6 threshold; most lab files are below it), every numeric value of that file becomes 0, SEQN included. The file is then cached for the rest of the process. The SEQN merge matches nobody, so every variable from that file is missing for that whole cycle, and complete-case analysis drops the cycle. No error is raised.

Direct test (`DEMO.xpt` 1999-2000, 9965 × 144 numeric cells, the same `mask` call repeated in one process): numexpr on, 2/60, 1/150 and 1/200 calls returned all zeros; numexpr off, 0/150. The comparison `abs() < 1e-70` alone was correct in 150/150; the failure is in the `where` step.

## 2 Evidence that this is what happened in the recorded M1 run

- The recorded N is the correct N minus exactly the 2007-2008 cycle in both papers: 9835a67923 3415 − 909 = 2506; 69700c2b09 3636 − 648 = 2988 (909 and 648 are the 2007-2008 complete cases in the correct sample).
- Both papers ran in the same process (batch 1) one minute apart (08:12 and 08:13), so one bad cached 2007-2008 file (lib cache of `DR1TOT_E`, which gives the dietary weight to both builders) explains both.
- A traced rerun of 69700c2b09, 9835a67923 and 71eedfeafb in one process on 15/11 (`tr_a3`, scratch) reproduced the recorded rows exactly: 69700c2b09 N 2988, 9835a67923 N 2506, both missing 2007-2008, while the file reads themselves were identical to the good runs (same rows, columns, missing count and sum). Two later traced runs caught the cached lib copy of `DEMO.xpt` 1999-2000 with sum 0 and no missing values (the good copy has 86,316 missing values and sum 2.916e10); both runs lost the 1999-2000 cycle of 71eedfeafb (N 28275 instead of 30410).
- Other one-cycle losses seen today with the unfixed code: 69700c2b09 lost 2011-2012 once (N 3050); a full single-process rerun lost 2005-2006 in four papers that read the dietary files through the same cache (71eedfeafb 27705, aca6f864fb 24533, 9835a67923 2840, a906145f51 5117).
- Ruled out: builder, library, spec table and data changes (source mtimes and the `.pyc` headers match the sources used at 08:12; all 337 downloaded files match the sha256 in `download_log.json`; nothing in `data-external/nhanes-extra/` changed after 07:56:43); `NHANES_ROOT` (not set); missing files (no existence check for an XPT file failed and then succeeded on recheck in any traced run); builder order (the same order gave both results).

Why the earlier replay "in one process, in the M1 batch order" gave the new N: the failure is random per large-file read (roughly 1 % to 3 % per call), not tied to the order.

## 3 Which N is correct

The larger N. Both papers describe samples that include 2007-2008: 9835a67923 "5 cycles of data from NHANES 2001–2010" (fulltext_txt/9835a67923.txt, line 22; reported N 3218), 69700c2b09 "NHANES spanning from 2005 to 2016" (line 13; reported N 4177). Neither builder excludes 2007-2008, and no paper in the sample drops a middle cycle. The recorded N lacked a cycle that the paper's own sample contains.

## 4 Fix

`tools/m1_lib_i139.py` now sets `pd.set_option("compute.use_numexpr", False)` at import. `m1_run_i139.py`, `m1_check_i139.py` and `m23_common_i139.py` (and through it M2, M3, M4, Jaccard, the any-combination and failure-cause scripts) import this module, so all later runs use plain numpy. The locked `nhanes_core.py` is not edited (KD3 hash). Without numexpr the good-run numbers are unchanged (section 6).

## 5 Old rows (as recorded at 08:12 to 08:13)

| rec_id | rank | OR reported | OR repro | 95 % CI repro | N reported | N repro | ΔN rel | Δ log OR |
|---|---|---|---|---|---|---|---|---|
| 9835a67923 | fail | 0.78 | 0.7646 | 0.5227-1.118 | 3218 | 2506 | 0.2212 | 0.02 |
| 69700c2b09 | fail | 1.07 | 1.066 | 0.9812-1.157 | 4177 | 2988 | 0.2846 | 0.004 |

## 6 New rows and check of every paper

Rerun with the fixed code at 19:0x by `python -I tools/m1_run_i139.py 9835a67923 69700c2b09` (log `extraction/m1_run_log_correction_20261115.txt`). The table was sha256 `cbeec4fb…53b1a706` and is now `3e8e5be0…ffd71b46`; only these two rows changed.

| rec_id | rank | OR reported | OR repro | 95 % CI repro | N reported | N repro | ΔN rel | Δ log OR |
|---|---|---|---|---|---|---|---|---|
| 9835a67923 | strict_pass | 0.78 | 0.7516 | 0.5517-1.024 | 3218 | 3415 | 0.0612 | 0.0371 |
| 69700c2b09 | near_pass | 1.07 | 1.046 | 0.9709-1.128 | 4177 | 3636 | 0.1295 | 0.0222 |

Other papers. Two full M1 reruns of all 47 papers with the fix, run at the same time as two processes (the setting in which the fault appeared before), each agree with the recorded table in 45/47 rows on status, rank, OR, CI, N, Δ log OR and ΔN; the two differing rows are the two above, with the values above. The other 40 testable papers (and the 5 not testable) keep their recorded values. Three of them (71eedfeafb, aca6f864fb, a906145f51) did lose a cycle in a rerun with the unfixed code today, but their recorded rows are the correct ones (confirmed by the fixed reruns), so the recorded M1 run was hit only in the two papers. The recorded value of the any-combination run (`m1_anycomb/`) and of the failure-cause run (`m1_fail_causes_i139.tsv`, `n_cc_check`) for the first combination equal the correct N in all 42 papers.

## 7 Effect on E1 (protocol §6, three denominators, Wilson 95 % CI)

| Quantity | Before | After |
|---|---|---|
| Strict pass | 9/42 = 0.2142 [0.117; 0.3593] | 10/42 = 0.238 [0.1348; 0.3852] |
| E1(i), testable | 15/42 = 0.3571 [0.2298; 0.5083] | 17/42 = 0.4047 [0.2704; 0.555] |
| E1(ii), not testable = fail | 15/47 = 0.3191 [0.2039; 0.4616] | 17/47 = 0.3617 [0.2396; 0.5046] |
| E1(iii), not testable = pass | 20/47 = 0.4255 [0.2951; 0.5671] | 22/47 = 0.468 [0.3333; 0.6076] |
| OR only, ≤ log 1.10 / ≤ log 1.25 | 29/42; 36/42 | 29/42; 36/42 (unchanged) |
| Pilot excluded, E1(i) / (ii) / (iii) | 14/40, 14/45, 19/45 | 16/40 = 0.4 [0.2634; 0.554]; 16/45 = 0.3555 [0.2321; 0.5016]; 21/45 = 0.4666 [0.3293; 0.6092] |
| Any combination passes | 18/42 | 18/42 (the any-combination run already had the correct samples) |

Reading under prereg §9.1 is unchanged: the upper bound of E1(i), 0.555, is below 0.80, so "there is a problem".

## 8 Downstream files changed by the correction

- `extraction/m2_results_i139.tsv`: M2 runs on M1 strict and near passes, so the two papers were added (two rows appended; the 15 earlier rows are kept as recorded). Both are "cannot be distinguished"; E2 width ratios 0.7324 (9835a67923) and 0.7511 (69700c2b09), both below 0.80. sha256 was `0c995293…31f7042`, now `173af2b1…4a2017`.
- E2: 3/14 = 0.2142 [0.0757; 0.4758] becomes **5/16 = 0.3125 [0.1416; 0.5559]**; median N-adjusted width ratio 0.9964 becomes 0.9839 (range unchanged, 0.6023 to 1.4965). Under prereg §9.1 the E2 reading moves from "inconclusive" to "there is a problem" (lower bound 0.1416 > 0.10). Pilot excluded: 2/13 becomes 4/15 = 0.2666 [0.1089; 0.5195].
- A full M2 rerun with the current code (not written to the record) also changes three earlier rows (6b7f67827b width ratio 1.03 to 0.9991, 94132d0c1a 1.0117 to 1.0013, a906145f51 0.9863 to 0.9992; the regime (i) n differs). These come from code edits after the M2 run (`tools/m23_common_i139.py` changed at 09:05, M2 ran 08:37 to 08:39), not from the numexpr fault (the n is larger, not smaller, in two of them). None crosses 0.80, so E2 is 5/16 either way. Not changed in the record; to be looked at separately.
- `extraction/m3_summary_i139.json` regenerated by `tools/m3_summary_i139.py` (no model refit): only `E3a_A_near_pass` and `E3a_B_near_pass` change, 0/15 [0; 0.2038] to 0/17 [0; 0.1843]. M3 itself already used the correct samples (n_paper_sample 3415 and 3636).
- `extraction/sens_sample_i139.json` regenerated by `tools/sens_sample_i139.py` (no refit): only the E1 and E2 entries above change.

## 9 Not checked

M2, M3, M3 sensitivity runs, M4 and the Jaccard sets were all produced with numexpr on and can carry the same one-cycle loss. Their n per paper does not show it directly (they use the correct weight, so n differs from M1 by design). A rerun of those steps with the fixed `m1_lib_i139.py` (imported by `m23_common_i139.py`) and a comparison of n per paper is needed before their numbers are final.

Update 15/11 22:0x: done; see `extraction/m234_correction_20261115.md`. The three M2 shifts named in section 8 come from the numexpr fault (one cached 2003-2004 MCQ file lost in the 08:37 M2 process), not from the 09:05 code edit.

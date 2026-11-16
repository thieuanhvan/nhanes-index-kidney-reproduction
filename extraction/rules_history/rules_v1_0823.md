# I139 M2 and M3: implementation rules (locked 15/11/2026 before any M2 or M3 model is run)

These rules make prereg_i139_v7 Section 6 (M2, M3) and Section 7 executable. Where the prereg is silent, the choice
below is made now, before results, and is reported as an implementation decision. Code: `tools/m3_run_i139.py`,
`tools/m2_run_i139.py`. Builders `tools/m1_specs/<rec_id>.py` give the paper's specification (first combination).

## Common to M2 and M3
1. **Correct weight** ("weight of the smallest subsample among the variables used"):
   - dietary index (rec_id in DIET): day-1 dietary weight WTDRD1 (2017-2020 WTDRD1PP; 1999-2002 pooled WTDR4YR),
     pooled with the NCHS factor (2 / total years; 4 for 1999-2002; 3.2 for 2017-2020);
   - else, if at least 99 % of the paper's complete cases have a positive fasting weight: pooled fasting weight;
   - else pooled MEC weight. Rows with weight 0 under the correct weight are dropped (subsample rule).
   Taylor linearisation, design df = #PSU - #strata, domain estimation (`nhanes_core.svy_logit`).
2. **Scale**: a continuous target is standardised by its weighted mean and SD in the analysed sample (log OR per SD);
   a 0/1 target (quantile contrast, cut point) is kept as the paper's contrast.
3. **Flag "index contains a component determined by the outcome"** (prereg Section 2): by code from
   `index_components`: serum albumin, haemoglobin, creatinine, cystatin C, BUN, urinary albumin, phosphorus, potassium,
   or uric acid when the outcome contains eGFR (CKD, DKD, low eGFR). Flagged papers stay in all analyses; the
   primary E3a denominator excludes them.

## M3 (Q1, main)
4. **Concept codings for variables added** (one fixed coding per concept, NHANES columns):
   age `age`; sex `female`; race RIDRETH1 5 groups (ref NH White); education 3 groups (ref < high school); PIR
   continuous; smoking never/former/current (ref never); diagnosed hypertension BPQ020 (1 yes, 2 no); diagnosed
   diabetes = DIQ010 1 or insulin DIQ050 1 or pills DIQ070/DID070 1 (0 when DIQ010 is 2/3 and no medication);
   physical activity = any leisure-time moderate or vigorous activity (GPAQ PAQ650 or PAQ665 = 1, 2007+; PAD200 or
   PAD320 = 1, 1999-2006; 0 when both answered no); alcohol = any drinking day in the past 12 months (ALQ120Q > 0 or
   ALQ121 1-10; skipped never-drinkers ALQ110/ALQ111 = 2 count 0); lipid-lowering drugs = BPQ100D 1 (0 when BPQ100D 2,
   or skipped because BPQ090D, BPQ080 or BPQ060 = 2); BMI/waist = BMI continuous. Codes 7/9 missing.
   If a concept's term is already among the paper's terms, it is not added twice (logged).
5. **Common sample**: domain of the builder, complete on y, the paper's terms, the variables missing under DAG A and
   under DAG B (`extraction/labels_i139.tsv`), correct weight > 0.
6. Models on the common sample: P (paper's terms), A (P + missing A), B (P + missing B). Δ_A = logOR(A) - logOR(P),
   Δ_B likewise. Sample-change difference = logOR(P on P's own complete cases) - logOR(P on common sample), reported
   separately.
7. **Paired Rao-Wu bootstrap**, B = 500, seed `20261114-i139-m3`: per stratum of the pooled design draw n_h - 1 PSUs
   with replacement, weight x n_h/(n_h - 1) x times drawn; the same draw for every model of a paper. Point fits by
   weighted IRLS. Check before: no stratum with fewer than 2 PSUs. CI = 2.5th and 97.5th percentiles; failed
   replicate fits are dropped and counted.
8. **Classification**: firm positive when the 95 % CI of Δ lies entirely outside [-log 1.10; log 1.10]; point-estimate
   positive when |Δ| > log 1.10 but the CI touches the interval; negative otherwise or when nothing is missing
   (no run needed, but controls are run for every testable paper).
9. **Negative control**: per paper a N(0,1) variable from seed `20261114-i139-negctl:<rec_id>`, set missing at random
   with the PIR missing rate of the builder domain; Δ_neg = logOR(P + rand) - logOR(P) on the common sample with rand
   observed. Firm positive proportion must be <= 0.05, else the threshold becomes the 95th percentile of |Δ_neg|
   (pre-specified deviation).
10. **Positive control**: Δ_pos = logOR(P) - logOR(P without the age terms) on the common sample; papers without an
   age term are skipped (counted).
11. **E3a**: stratified weighted proportion of firm positives with stratum weights N_h* = N_h x full-text eligibility
   proportion of stratum h; variance with finite population correction, logit CI; stratified bootstrap by paper
   (2000 draws, percentile) when a stratum is 0 % or 100 %. Also unweighted proportion with Wilson CI.
   Primary: DAG A, testable and unflagged. Parallel: DAG B; secondary denominators per prereg Section 7.
12. Descriptive analyses not in this run (later, not used for the main conclusion): E3b, two-order decomposition,
   PR/marginal sensitivity, other sensitivity analyses, RoE.

## M2 (Q3, near-pass papers: strict + near, 15)
13. Paper's specification on the builder domain. Regimes: (i) correct weight, design SE; (ii) wrong-type weight
   (correct fasting or dietary -> MEC; correct MEC -> interview), design SE; (iii) unweighted, model SE;
   (iv) correct weight, model SE; (v) unweighted, design SE. Model SE = inverse information with weights normalised
   to mean 1, normal 95 % CI.
14. For each regime: |Δ log OR| to the reported OR and |Δ log CI width| with the reproduced width (log hi - log lo)
   multiplied by sqrt(n_reproduced / n_reported). Matching regime = closest on both components and at least log 1.10
   from the next regime on at least one component; else "cannot be distinguished".
15. E2: width ratio = reported log width / (regime (i) log width x sqrt(n_reproduced / n_reported)); proportion < 0.80
   among the 15, with median and per-paper values. Effect size: ratio of ORs per SD (matched regime vs (i)).

## Amendment 1 (15/11/2026, after the code test on one paper, f1ad3fb0ff; no other M3 result seen)
Rule 1 threshold "at least 99 %" changed to "at least 80 %". Reason (data structure, not results): NCHS gives a zero
fasting weight to about 5 % of people who have fasting-session lab values but did not meet the fasting criterion, so
papers that use fasting glucose/TG reach only 0.90-0.97, while papers without fasting variables sit at 0.43-0.52
(share of complete cases with positive fasting weight, all 42 testable papers). 80 % separates the two groups.
The f1ad3fb0ff test output (MEC weight) is discarded and rerun.

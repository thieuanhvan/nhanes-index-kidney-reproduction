# Full pre-registration I139, version 7 (14 November 2026)

Thieu Anh Van, University of Information Technology, VNU-HCM

Working title: Confounder specification in cross-sectional NHANES studies of an index, ratio or score and kidney disease: a per-paper reanalysis of a stratified random sample drawn from a Europe PMC frame.

Version 7 (14 November 2026). English translation of version 6 (`prereg_i139_v6_vi.md`, sha256 `72db630ce54c3d65d90c8d175aaf0eac10145d967020fc802387d147fd25ca21`), the version posted on OSF, plus three pre-specified edits that do not change any outcome:

- **Edit A:** a short "Related work and novelty" subsection after Section 1, with a narrowed novelty claim.
- **Edit B:** Appendix A now records the sha256 of the DAG folder.
- **Edit C:** the 9-item NCHS reporting checklist in Section 7 is replaced by a 12-item list reconciled with the NCHS guidelines of 26 January 2026.

Versions 1 to 6 and the pilot results are kept unchanged for comparison. Version 6 revised version 5 (`prereg_i139_v5_vi.md`, sha256 `20bdfbf6a952f03d5b791a6d5de7e16cd2c2409361289ac7c372dbc3dc34207b`) in response to the blinded round-2 review `reviews/prereg_i139_v5_review_vi.md` (6/10, MAJOR REVISION); the accept / partly accept / reject table for each point is in Section 13. The largest change in version 6: the main question moved entirely to **confounder specification** (omitted confounding, mediator adjustment); E3a is the only primary outcome; reproduction (E1) and survey design (E2) became key secondary outcomes. The sha256 hash of this file is stored in the `.sha256` file of the same name. This file is posted on OSF before any paper in the sample is extracted.

**Scope of inference.** Every proportion in this document generalises to "the eligible papers in the frame of 282 Europe PMC records that mention NHANES, an index/ratio/score and kidney disease in the title", not to the whole NHANES literature. The target estimate of E3a is **frame-wide** (all index families). Results for TyG alone (stratum F1, expected 9 to 12 papers) are reported only as counts x/y with a per-paper description, not as proportions with confidence intervals.

## 0. The pilot and what it must not be used for

The pilot (7 TyG × kidney disease papers, `pilot/`) was used only for design. After fixing four technical errors (the prereg name in the source code, the flag for same-as-paper specification, missing variables coded as 0 or as a category, serum creatinine of 2005 to 2006 not yet calibrated), the main pilot numbers did not change: of a grid of 189 specifications, 82.53% were statistically significant; 54.76% (CKD, adults) and 71.43% (albuminuria, adults) are proportions within the two covariate sets that contain diabetes, hypertension and BMI; the unweighted version gave 98.94%; reproduction succeeded for 3/4 papers.

An exploratory analysis done after the fixes (not pre-registered): with the minimal covariate set (age, sex, race/ethnicity, education, PIR, smoking), 63/63 specifications were statistically significant, with ORs per 1 SD of TyG from 1.312 to 1.537. The "fragility" seen in the pilot came mainly from over-adjustment, not from errors in the literature. This is why the main question moved to decomposing the effect size by covariate set, and why E3b (mediator adjustment) is reported descriptively only and is not a primary outcome.

Pilot papers B1 to B7 are inside the frame (stratum F1). If selected, they are extracted and reanalysed from scratch with the procedure of this version. Sensitivity analysis: every outcome is recomputed with the pilot papers excluded.

Choices made after seeing the pilot data (declared so that readers can judge for themselves): dropping the cycle-window grid; the DAG-based minimal set; opening the k grid to percentiles 1 to 99; the log 1.10 and log 1.25 thresholds (present since version 3). The pilot numbers quoted here are design information that came from data already seen, not independent evidence.

## 1. Research questions and primary estimand

Within the stated frame:

- **Q1 (main question, confounder specification):** do the published papers omit certain confounders (according to the pre-declared DAG) to the point that adding those confounders changes the effect size by at least 10%? Do the papers adjust for components or mediators of the index itself, and by how much does this shift the OR?
- **Q2 (computational reproduction):** from the public CDC files and the description in the paper, can the main estimate be reproduced?
- **Q3 (survey design):** does the paper handle subsample weights, strata and PSUs correctly, and when corrected, how much do the effect size and the width of the confidence interval change?
- **Q4 (inflection point):** for papers that report a threshold or inflection point, can the data locate that threshold, given a method shown to locate a true threshold and not to locate a false one when none exists?

Every DAG reference remains a cross-sectional association, not a causal effect.

Effect-scale convention: the OR per 1 SD is computed on exactly the scale the paper used (SD of ln index if the paper used ln; weighted SD within the paper's analytic domain); the OR per paper unit is reported alongside. For papers that report categories: every comparison uses the paper's own cut points (for example Q4 versus Q1); the OR per SD is supplementary only.

### 1.1 Related work and novelty (Edit A)

Related work found in the novelty search of 14 November 2026 (Europe PMC, CrossRef, arXiv; Semantic Scholar returned HTTP 429 once and was dropped; the search log is in `search_log_vi.md` and the gap matrix in `gap_matrix_i139_vi.md`):

- [Ahern et al. 2025, Epidemiology](https://doi.org/10.1097/ede.0000000000001821) changed the adjustment set of 33 BMI-mortality estimates to a sufficient minimal set, using NHANES III as external data and bias factors; the pooled result barely moved (0.88 to 0.90) and heterogeneity fell.
- [Luo et al. 2026, Int J Hyg Environ Health](https://doi.org/10.1016/j.ijheh.2026.114758) examined 175 air-pollution cohorts and found that 69.7% may be over-adjusted, with BMI the most frequent mediator.
- [Mekonnen et al. 2026, Crit Rev Food Sci Nutr](https://doi.org/10.1080/10408398.2026.2719168) reviewed 46 ultra-processed food papers: 2 used a DAG and 38 adjusted for mediators.
- [van Zwieten et al. 2024, IJE](https://doi.org/10.1093/ije/dyad177) studied how systematic reviews handle over-adjustment.
- [Patel, Burford, Ioannidis 2015, J Clin Epidemiol](https://doi.org/10.1016/j.jclinepi.2015.05.029) introduced the ratio of extremes (RoE) for the vibration of effects under different adjustment sets; it is used here for the RoE sensitivity analysis (Section 6, M3).

Narrowed novelty claim: this study is a per-paper reanalysis of published NHANES studies on their own public data, which separates omitted confounding (E3a) from mediator or component adjustment (E3b), for index × kidney outcomes. Counting over-adjustment is already well covered by the work above, which is why E3b is reported descriptively only. No study was found that does this per-paper reanalysis on at least 10 NHANES index × kidney papers, so the KD5 trigger (b) in Section 9.2 was not met on 14 November 2026.

## 2. Frame, eligibility criteria and sampling

### 2.1 Frame (locked)
- Source: Europe PMC, run on 14 November 2026, 362 records. Query P3, verbatim:
  `(TITLE:albuminuria OR TITLE:"albumin-to-creatinine" OR TITLE:"albumin creatinine" OR TITLE:"diabetic kidney disease" OR TITLE:"diabetic nephropathy" OR TITLE:"chronic kidney disease") AND (TITLE:index OR TITLE:ratio OR TITLE:score) AND (NHANES OR "National Health and Nutrition Examination Survey")`
- Keep records that contain "NHANES" or "National Health and Nutrition Examination Survey" in the title or abstract; remove duplicates by DOI, then by title: **282 records**.
- Frame file: `frame/frame_i139.tsv`, sha256 `b5d7db0e4c9639946b48a628d56e51d193b21fb05f17b6b81ecc8c5a31a55872`, built with `tools/build_frame_i139.py` (rerunning gives the same hash).
- Strata (assigned by title keywords in a fixed order in the source code): F1 TyG and derivatives (15), F2 inflammation and haematology (53), F3 diet and lifestyle (79), F4 lipids, metabolic, anthropometric (55), F5 other (80). The stratum assignment only affects allocation; the covariate roles (Section 5) follow the **DAG family** of the index, assigned by code from the list of extracted components (Section 5.3).
- The frame is not extended after registration.

### 2.2 Eligibility criteria
A paper is eligible when all five of the following hold:
1. It uses NHANES data (cycles 1999 to 2023) as the main data, and the NHANES part can be separated if the paper has several sources.
2. The main exposure is an index, ratio or score computable from public NHANES variables.
3. The outcome is prevalent kidney disease: CKD, albuminuria, low eGFR, or CKD in people with diabetes (called DKD by the paper).
4. There is at least one adjusted association estimate (OR, PR or β) with a confidence interval or standard error.
5. The analysis is cross-sectional (or has a separable cross-sectional part).

Excluded: mortality or progression outcomes in people who already have CKD; UACR or eGFR as the exposure; comments, reviews, pooled analyses; retracted papers; indices that need restricted NCHS data.

**Flag "index contains a component determined by the outcome"** (set at screening, the paper is not excluded): the index contains creatinine, cystatin C, BUN, urinary albumin, uric acid (when the outcome contains eGFR), haemoglobin, serum albumin, phosphorus or potassium. Every outcome is reported separately for the flagged and unflagged groups; the primary estimate of E3a uses only the unflagged group.

**Descriptive label "conditioning on a descendant of the exposure":** attached to a paper that restricts the sample to people with diabetes (DKD) when the index contains glucose, insulin or HbA1c (families D1, D4b in Section 5). It does not change the analysis; the number of papers is reported.

### 2.3 Sampling
- Seed: the string `20261114-i139-v4` (kept from version 4 because the frame and ordering were already locked). The random order within each stratum is the order of sha256(seed + ":" + record key), recorded in the `rank_in_stratum` column of the frame file.
- Step 1: screen all 282 records for eligibility at title and abstract level (Section 3), giving the number of eligible papers N_h in each stratum.
- Step 2: within each stratum, follow `rank_in_stratum`, take papers eligible at abstract level, and read the full text. A paper that is not eligible at full-text level is logged with the reason and the next paper is taken. A paper without full text: try through the university library and by emailing the corresponding author within 7 days; only after that is it replaced by the next paper.
- Allocation (target 50 papers, acceptable 40 to 60): F1 takes all eligible papers (expected 9 to 12); F2, F3, F4, F5 take 10 papers each. If a stratum falls short, take all of it and shift the shortfall to other strata in the order F4, F2, F5, F3, one paper at a time. Weights are always computed from the actual n_h.
- Extrapolation weight of a paper in stratum h: N_h* / n_h, with N_h* = N_h × (the eligibility proportion at full-text level observed in stratum h).
- Report by stratum the number of papers replaced for lack of full text. Sensitivity analysis: treat these papers as "not testable" instead of replacing them.
- Expected precision (illustrative example, eligibility 50% in every stratum, so F1 has about 7.5 papers, the other strata 10 papers, n = 47.5): weight 1 in F1 and 2.650 to 4.000 in the other strata give a Kish factor of about 1.121, an effective sample size of about 42.35, and a 95% CI half-width around 50% of about ±0.1505 when the finite population correction is ignored and about ±0.1259 with it (F1 is taken in full, so its sampling fraction is 1). With 30 papers, the unweighted half-width is already ±0.1789. Consequence for E3a: if the denominator is only about 25 papers, 0/25 still gives a Wilson upper bound of 0.1331; about 40 papers with 0 positives (upper bound 0.08762) are needed to bring the upper bound below 10%. A "negative" conclusion for E3a is therefore only "no evidence seen, precision limited" (Section 9).

## 3. Screening, extraction and verification (declaration of AI tool use)

### 3.1 Declaration
Screening, extraction and reconciliation are carried out by agents based on large language models (LLMs), under the supervision of the single author. There is no second human rater. The author personally checks the three fields that decide E3a following step 10 below; this is reported as an "author audit", not as inter-rater reliability. Agreement between machine passes is reported as "agreement between two LLM extraction passes". The models used are Claude (Anthropic; Opus and Sonnet models, via command-line and web interfaces); model names, versions, run dates and the verbatim prompts are kept in the log and released with the data. The reproduction part is deterministic source code, with no manual judgement beyond the locked specification table.

### 3.2 Verification procedure
1. **Two independent passes.** Screening and extraction are each run twice in two separate sessions, neither pass seeing the result of the other, with two independently written prompts. For the fields that decide E3a (the model of the main estimate, the covariate list of that model, the OR and CI), the two passes **use two different model families**. If access to a second family is not available, this is stated in the paper, and agreement between the two passes is not used as evidence of accuracy (then only step 5 and step 10 are used). Report Cohen's κ for screening; for extraction, report per field the raw agreement, Cohen's κ and Gwet's AC1 (less inflated when one value dominates).
2. **Verbatim evidence.** Every key field (cycles, domain, N, outcome definition, eGFR equation, UACR threshold, covariate set, weight variable, how cycles were pooled, main estimate and CI, which model is the model of the main estimate, components of the index, inflection point if any) must carry a verbatim excerpt and its location (table/paragraph). The covariate set must be taken both from the table footnote and from the methods section; if the two differ, record both and use the footnote of the table that holds the main estimate.
3. **Automatic checks** (source code, run on every extraction):
   - the excerpt is present verbatim in the full text (string match after normalising whitespace and dashes);
   - the covariate set and the main estimate are taken from the same table and the same model (location check);
   - numeric consistency: the OR lies inside the CI; the CI is symmetric on the log scale within rounding error; p (if reported) matches the CI;
   - N does not exceed the number of people with data in the CDC file for the exact cycles, age group and subsample the paper describes;
   - variable names, weight variable and cycle exist in the column catalogue of the corresponding XPT file;
   - the weight variable the paper names is the right type for the smallest subsample among the variables it uses.
4. **Reconciliation.** A field with disagreement, or one that fails a step-3 check, goes to a third pass: read the full text, quote the verbatim passage, and record the reason in the log. Cases where both passes make the same error do not enter reconciliation; that is the reason for steps 5 and 10.
5. **Controlled error injection (60 errors).** Take the reconciled version of 15 papers (chosen with the seed `20261114-i139-inject`), and plant by code 60 known errors, each **keeping a real excerpt from the full text** (no error is caught merely because the excerpt does not match):
   - 30 "crude" errors: wrong OR or CI, wrong N, wrong cycle, wrong weight variable, with a real excerpt that does not contain the planted value;
   - 30 "plausible but wrong" errors, with a real excerpt that matches the planted value: switching to another model in the same table (for example Model 2 instead of Model 3), to a sensitivity model, to the covariate list of another model in the same paper, dropping a covariate mentioned only in the methods section, or to the OR of a different outcome.

   Rerun step 3 and an independent extraction pass (a different model family from the pass that produced the reconciled version, if available). Report the overall detection rate and the rate in the two groups, with Wilson CIs. Condition to proceed: the lower bound of the 95% CI of the overall detection rate is at least 0.85 (equivalent to at least 57/60, lower bound 0.8629) and the point estimate in the "plausible but wrong" group is at least 0.85 (at least 26/30).
6. **Prompt-tuning set separate from the test set.** The prompt is tuned on 4 pilot papers (chosen by seed), then locked. The test set consists of the 3 remaining pilot papers plus 5 random papers taken from the tail of `rank_in_stratum` that is not selected into the sample. A **reserve test set** of the next 5 papers in that tail is locked in advance, used only for the retest after a prompt revision (the first test set is not reused). Condition to proceed: agreement of at least 0.90 between the reconciled version and the full text (checked through steps 2 and 3) on the key fields, and at most 5% of fields failing the verbatim check.
7. **The specification table and the "first combination" (M1)** are produced by the extraction pass before any model of that paper is run; the log is time-stamped and the table is hashed with sha256 before the first run.
8. **N_h from LLM screening.** Report the screening κ; sensitivity of the outcomes when N_h is computed from each screening pass separately rather than from the reconciled version.
9. **Variance implementation check.** Before M1, the survey-estimation source code must give the same coefficients and standard errors (difference at most 1%) as an independent implementation (R package `survey`, or `samplics` if R is not available) on 3 pilot specifications, and must reproduce an estimate given in the CDC NHANES tutorial.
10. **Author audit (mandatory, narrow).** After the specification table is locked and the covariate labels are assigned (Section 5.4), **before any M3 model is run**, the author checks against the full text exactly three fields: which model is the model of the main estimate; the covariate list of that model; the OR and CI. Scope: (a) every paper with the label "missing at least 1 confounder" under DAG A or DAG B; (b) 10 papers chosen with the seed `20261114-i139-audit` from the remaining papers. The audit result is hashed with sha256 before the first M3 run. Report accuracy per field with Wilson CIs. Condition to proceed (into KD2): for the covariate list field, the number of papers in which an extraction error changes the "missing confounder" label is at most 10% of the audited papers; if exceeded, this field is re-extracted for every paper in the sample under the author's checking, and only then is M3 run. For papers that "nearly pass" in M1, reproducing the reported estimate is additional objective evidence that the model and covariates were extracted correctly; this proportion is reported separately.

## 4. NHANES data and index reconstruction
- XPT files are downloaded directly from CDC when needed, and the sha256 of each file is written to the log; files processed by third parties are not used.
- Serum creatinine is calibrated following the NHANES recommendation for cycles 1999 to 2000 and 2005 to 2006 before eGFR is computed.
- Missing values (including codes 7/9 "refused/don't know") are left missing; complete-case analysis as in the papers; no "missing" category is created.
- Indices are reconstructed with the formula the paper gives. Indices that need a food composition database (for example HEI, DII) are reconstructed from the public dietary files if the paper names all components; otherwise the paper is classed "not testable (index cannot be reconstructed)".
- Cycles 2017 to March 2020 (pre-pandemic) and 2021 to 2023 use their own weights following NCHS guidance.
- Rao-Wu bootstrap: the public NHANES files have 2 PSUs per stratum, so each draw takes 1 PSU with replacement in each stratum and rescales the weights following Rao-Wu; check beforehand that after pooling cycles no stratum is left with only 1 PSU.
- Overlap of people between papers: Jaccard index on the set of SEQN of the reconstructed analytic sample, computed for each pair of papers with overlapping cycles.

## 5. DAG, node table and covariate labels

### 5.1 Two pre-declared DAGs
- **DAG A (main):** the minimal set contains only the nodes ranked "definite confounder" (G) in table 5.2.
- **DAG B (alternative, reported in parallel):** DAG A plus the nodes marked "B+" in table 5.2 (BMI/waist circumference as a confounder when it is not a component; alcohol; physical activity; lipid-lowering drugs in families containing TG; diagnosed diabetes in the pure-lipid family).

DAG A is conservative (few nodes, E3a hard to make positive); DAG B is close to the usual way of adjusting in the literature. E3a is computed under both (Section 7). The DAG figures and dagitty code for each family are in Appendix A of this version, written from table 5.2 (the table is the only source; the dagitty code is generated from the table by source code).

### 5.2 Node table by DAG family

Symbols: **G** definite confounder (in DAG A and B); **T** contested (not in DAG A; in DAG B if marked **B+**; otherwise used only in sensitivity analysis); **K-C** do not adjust because it is a component of the index; **K-M** do not adjust because it is a mediator or descendant of the exposure; **K-O** do not adjust because it is a descendant or component of the outcome (conditioning on the outcome); **K-X** do not adjust because it is a collider; **P** precision variable, included if the paper has it but not counted toward the "missing" label.

DAG families: **D1** TyG and derivatives (TyG, TyG-BMI, TyG-WC, TyG-WHtR); **D2** inflammation and haematology; **D3** diet and lifestyle; **D4a** pure lipid (for example AIP, TG/HDL, non-HDL/HDL); **D4b** indices with glucose, insulin or anthropometry (for example HOMA-IR, METS-IR, VAI, LAP, WHtR); **D5** other (environment, sleep, mental health, social indices, etc.).

| Node | D1 | D2 | D3 | D4a | D4b | D5 | One-line reason |
|---|---|---|---|---|---|---|---|
| Age, sex | G | G | G | G | G | G | Affects both the index and kidney function, not a descendant |
| Race/ethnicity | G | G | G | G | G | G | Black people have lower TG and higher CKD risk (negative confounding) |
| Education, PIR | G | G | G | G | G | G | Social position affects diet, metabolism, access to treatment and kidney disease |
| Smoking | G | G | G | G | G | G | Affects lipids, insulin resistance, inflammation and kidney damage; for D3 with smoking inside the index it is K-C |
| Alcohol | T, B+ | T, B+ | T, B+ (K-C if in the index) | T, B+ | T, B+ | T, B+ | Affects TG, inflammation and the kidney; measurement changes between cycles, so not in DAG A |
| Physical activity | T, B+ | T, B+ | T, B+ (K-C if in the index) | T, B+ | T, B+ | T, B+ | As above |
| BMI, waist circumference, WHtR | T, B+ (K-C in TyG-BMI, TyG-WC, TyG-WHtR) | T, B+ | K-M | T, B+ | K-C or K-M | T, B+ | Obesity causes insulin resistance and glomerular hyperfiltration (confounder) but may also lie on the path from the index (mediator); kept in the contested rank. For D3, diet leads to obesity, so it is a mediator |
| Diagnosed diabetes (self-report, medication) | K-M | G | G | T, B+ | K-M | T | D1, D4b: diabetes is a descendant of the glucose component; D4a: diabetes causes dyslipidaemia and kidney disease but is also a descendant of insulin resistance; D2, D3: diagnosed people change diet, have chronic inflammation and high kidney risk |
| HbA1c, fasting glucose, insulin | K-C or K-M | K-M (use diagnosed diabetes instead of the measurement) | K-M (as D2) | T | K-C or K-M | T | Glycaemic measurements are descendants or components in families with glucose |
| Diagnosed hypertension | K-M | G | G | K-M | K-M | T | D1, D4: insulin resistance leads to hypertension; D2, D3: hypertension causes chronic inflammation, lifestyle change on advice, and kidney disease |
| Measured blood pressure | K-O | K-O | K-O | K-O | K-O | K-O | Cross-sectional: measured blood pressure is partly a descendant of CKD |
| Lipid-lowering drugs (statin, fibrate) | T, B+ | T | T | T, B+ | T | T | In families containing TG, drugs change the exposure value and are prescribed by cardiovascular risk related to CKD |
| Diabetes drugs (insulin, metformin, SGLT2i, GLP-1RA) | T | T | T | T | T | T | Lower glucose and TG; SGLT2i and GLP-1RA lower albuminuria; both a confounder and a descendant of diabetes. Separate sensitivity analysis |
| Antihypertensive drugs (ACEi/ARB, diuretics) | T | T | T | T | T | T | ACEi/ARB lower albuminuria; diuretics raise TG and glucose. Separate sensitivity analysis |
| HDL, LDL, total cholesterol | T (HDL recorded as a plausible descendant of insulin resistance) | T | T | K-C if a component, T if not | T | T | Not enough consensus to enter DAG A or B |
| CRP, white blood cells | T | K-C if a component, T if not | T | T | T | T | Inflammation may be a mediator or a confounder |
| Uric acid | K-O | K-O | K-O | K-O | K-O | K-O | Affected by kidney function |
| eGFR, creatinine, cystatin C, UACR, BUN | K-O | K-O | K-O | K-O | K-O | K-O | Measure kidney function or damage: the same process as the outcome, even when the outcome is only albuminuria or only eGFR |
| Diagnosed cardiovascular disease | K-X | K-X | K-X | K-X | K-X | K-X | Descendant of both the metabolic index and CKD |
| NHANES cycle | P | P | P | P | P | P | Time trend, changes in laboratory methods; not counted toward the "missing" label |

**Reverse causation and timing.** A cross-sectional design cannot establish temporal order. Nodes that may be descendants of CKD are marked K-O (measured blood pressure, uric acid). Diagnosed hypertension, CRP and HDL may also be affected in reverse by CKD; this is one reason they sit in K-M or T, not in G.

### 5.3 Assigning the DAG family
The DAG family is assigned by source code from the list of extracted components, in this order: contains TG and glucose, assign D1; contains glucose, insulin, HbA1c or anthropometric measures, assign D4b; contains only lipids, assign D4a; contains inflammatory or haematological markers, assign D2; contains dietary or behavioural components, assign D3; otherwise assign D5. The assignment is made and hashed before M3 is run. The stratum assignment F1 to F5 of the frame does not change.

### 5.4 Variable equivalence table and labels
Each variable a paper adjusted for is assigned by source code to a concept code following a pre-locked equivalence table. A concept counts as "present" if the paper has any of the forms below:

| Concept | Counts as present | Cannot substitute |
|---|---|---|
| Age | continuous or in groups | |
| Race/ethnicity | NHANES variable with at least 3 groups (RIDRETH1/RIDRETH3); binary "white/other" counts as present but is labelled "crude" (the sensitivity analysis treats it as missing) | |
| Education | any ordering with at least 2 groups | |
| PIR | PIR continuous or in groups; household or family income in bands; binary below the poverty threshold | health insurance, employment, marital status |
| Smoking | self-report (current/former/never, pack-years) or serum cotinine | |
| Alcohol | any measure of amount or frequency, including binary | |
| Physical activity | any measure (MET, minutes per week, categories) | |
| BMI/waist | BMI, waist circumference or WHtR | |
| Lipid-lowering drugs | statin use or self-reported "taking lipid-lowering medication", or from the medication file | self-reported high cholesterol |
| Diagnosed diabetes | self-reported physician diagnosis or diabetes medication | diabetes defined with glucose/HbA1c (separate code "diabetes by measurement") |

Each paper is labelled by source code: "sufficient for confounding (A)" / "missing confounder (A)", "sufficient (B)" / "missing (B)", with the list of missing variables; "adjusts for a component"; "adjusts for a mediator"; "conditions on the outcome"; "conditions on a collider". There is no per-paper manual judgement; cases the equivalence table does not cover are logged as deviations with the decision.

## 6. Analyses

Correct design in every analysis: the weight of the smallest subsample among the variables used, divided by the number of pooled cycles (or the NCHS formula for special cycles); Taylor linearisation variance with strata and PSUs; degrees of freedom = number of PSUs minus number of strata; t-based CIs; subpopulations by domain estimation. Sensitivity for degrees of freedom: df = design df + 1 − number of parameters (the `summary.svyglm` convention).

### M1. Specification reconstruction and reproduction (Q2)
- Main estimate of the paper: the most fully adjusted estimate in the abstract, preferring the continuous index. If the abstract has several outcomes: CKD first, then albuminuria, then low eGFR; within the same type, the first one stated.
- Ambiguity: at most 2 readings per ambiguous point, no more than 8 combinations in total; the "first combination" is the most literal reading, recorded in advance (Section 3.2, step 7).
- Ranking: **strict pass** |Δ log OR| ≤ log 1.10 and N difference ≤ 10%; **near pass** |Δ log OR| ≤ log 1.25 and N difference ≤ 15%; **fail**; **not testable** (reason recorded). Also report the cumulative curve of |Δ log OR|.
- Rounding: if the reported OR has fewer than 3 significant digits after subtracting 1, compare on the per-SD scale; or count as a pass when the reproduced estimate, rounded to the same number of digits, gives exactly the reported figure.
- Outcomes use the first combination. Sensitivity: "any combination passes", with the number of passing combinations over the number of combinations run.

### M2. Survey design (Q3), run only on "near pass" papers
- With the paper's specification, run: (i) correct design; (ii) wrong-type weights, design standard errors; (iii) unweighted, model standard errors; plus two split versions: weighted with model standard errors, and unweighted with design standard errors.
- **Regime the paper used:** report |Δ log OR| and |Δ log CI width| separately (width adjusted by the square root of the N ratio). A regime counts as matching only if it is the closest on both components and is at least log 1.10 away from the next regime on at least one component; otherwise it is classed "cannot be distinguished".
- Effect size: the ratio of ORs per SD (regime used versus correct design) and the ratio of CI widths on the log scale, with the width multiplied by the square root of (reconstructed N / reported N).

### M3. Confounder specification (Q1), main analysis
Run on every "testable" paper (specification reconstructed under the first combination), with correct design, on the log OR per SD scale (or per the paper's cut point).

- **Common sample.** For each paper, determine the set of people with complete data for the union of: the paper's set, the variables missing under DAG A, the variables missing under DAG B. Every model compared in M3 is run on this set. The difference due to changing the sample (the paper's set on the paper's sample versus the paper's set on the common sample) is reported separately and does not count toward E3a.
- **E3a comparison (separate, no decomposition chain):** Δ_A = log OR(paper's set + variables missing under DAG A) − log OR(paper's set), on the common sample; likewise Δ_B for DAG B. No variable of the paper is dropped in this comparison.
- **Uncertainty:** the 95% CI of Δ by paired Rao-Wu bootstrap (the same draw for both models), B = 500.
- **Classification of each paper:** **firm positive** when the 95% CI of Δ lies entirely outside the interval [−log 1.10; log 1.10]; **point-estimate positive** when |Δ| > log 1.10 but the CI touches that interval; **negative** when |Δ| ≤ log 1.10 or the paper is missing no variable. A paper missing no variable counts as negative without running.
- **E3b comparison (descriptive):** drop the K-C, K-M, K-O, K-X variables from the paper's set, on the common sample. Pre-declared expected direction: when these variables are present, the OR shrinks toward 1. E3b describes the paper targeting a different estimand (an "independent" association); it must not be called an error.
- **Two-order decomposition:** from the paper's set to the DAG A set by two routes: (add the missing variables, then drop the non-adjusted variables) and (drop first, add after). Report both, with the mean of the two routes. The full chain from the specification as reported, through reproduction and correct design, to the DAG A set is kept as in version 5, for description only.
- **Non-collapsibility of the OR:** the sensitivity analysis uses the prevalence ratio (PR) from a survey-weighted Poisson model and the marginal standardised estimate; report the number of papers that change classification.
- **Negative control for E3a:** add to the paper's set a random variable generated with the seed `20261114-i139-negctl`, with a missingness rate equal to that of PIR; the proportion of "firm positive" papers must be at most 0.05. If exceeded, the threshold is tightened to the 95th percentile of |Δ| under the negative control and recorded as a pre-specified deviation.
- **Positive control for E3a:** drop age from the paper's set and compute Δ when adding it back; report the proportion of "firm positive" papers to show that the comparison can detect a real confounder.
- **Other sensitivity analyses (descriptive only):** CKD-EPI 2009 versus 2021; UACR > 30 versus ≥ 30 mg/g; adding diabetes drugs or antihypertensive drugs; binary race treated as missing; the "conventional" set (DAG A plus diagnosed diabetes and hypertension) for every family; for albuminuria, a single UACR sample versus albuminuria confirmed by a second sample in cycles that have a second-sample file. Report the RoE (the ratio of the 97.5th to the 2.5th percentile of the OR per SD, following Patel 2015).

### M4. Inflection point (Q4), only for papers that report a threshold or inflection point
- (a) Estimate k: two-segment logistic model on the paper's specification, correct design; a grid of 100 points from the 1st to the 99th percentile with a weighted profile likelihood; report the proportion of bootstrap samples that pile up at the edge of the grid.
- (b) 95% CI of k: Rao-Wu bootstrap, B = 500; report in index units and in SD.
- (c) Nonlinearity test: design-based Wald on the coefficients of a restricted cubic spline (4 knots at percentiles 5, 35, 65, 95).
- (d) **Positive control:** simulate 200 outcome sets from a two-segment model with true k equal to the paper's k, at two slope-change levels (as in the paper, and 50% of that). Keep covariates, weights, strata-PSU and N; add a PSU-level random effect on the logit scale. Localisation power = the proportion of sets whose CI of k is narrower than 1 SD and contains the true k; coverage = the proportion of CIs containing the true k.
- (e) **Negative control:** simulate 200 sets from the fitted smooth spline model. False localisation rate = the proportion of sets whose CI of k is narrower than 1 SD.
- Budget: in simulation, k is estimated by an iterative Muggeo-type algorithm, with bootstrap B = 200 per set; real data keep the grid and B = 500. At most 2 parallel processes.
- Classification (valid only when coverage in (d) is at least 0.90 and false localisation in (e) is at most 0.10): **localisable, matching**; **localisable, different**; **not localisable** (CI at least 1 SD wide and power at the 50% slope-change level at least 0.80); **inconclusive**. Sensitivity: threshold of 0.5 SD.

## 7. Outcomes and estimation

**The single primary outcome: E3a(A)** = the weighted proportion of "firm positive" papers (M3, DAG A), with the denominator being the "testable" papers without the flag "component determined by the outcome". This is the only estimate used for the main conclusion.

Reported in parallel, with the reading pre-declared:
- **E3a(B):** as above under DAG B. If E3a(A) and E3a(B) give the same conclusion category (Section 9), the conclusion is said to be independent of the DAG; if they differ, the conclusion is recorded as "depends on the DAG choice" and both are stated. Report the number of papers that change label between A and B.
- Secondary denominators for E3a: (i) papers that "nearly pass" in M1; (ii) the whole sample, counting "not testable" papers with a "missing" label as negative (lower bound) and as positive (upper bound); (iii) the flagged group.
- **E3a-label** (descriptive, needs no NHANES data): the proportion of papers with the label "missing confounder" under A and under B, with a list of the commonly missing variables.
- E3a by point estimate ("firm positive" + "point-estimate positive") is a sensitivity analysis.

Key secondary outcomes:
- **E1:** the proportion of "near pass" papers (M1, first combination), with three denominators: (i) testable papers; (ii) the whole sample, counting "not testable" as fail; (iii) the whole sample, counting it as pass.
- **E2:** among "near pass" papers, the proportion with a CI width ratio (adjusted by N, M2) < 0.80, with the median and the per-paper distribution. Known limitation: the 0.80 threshold only catches a missing design when the design factor is about 1.562 or higher.

Descriptive outcomes: **E3b** (the proportion of papers that adjust for a component or mediator and the size of the OR shift), **E4** as a count x/y (papers "not localisable" among threshold papers with a valid classification), proportions by each adjustment label, the proportion "cannot be distinguished" in M2, the median RoE, overlap of people (Jaccard), the label "conditioning on a descendant of the exposure", the TyG result (D1) as x/y, and a reporting checklist following the NCHS NHANES analysis guidelines, rated yes/no/not applicable (Edit C, 12 items):
1. States the survey cycles (a STROBE reporting item, not in the NCHS text).
2. States the name of the weight variable.
3. Uses the weight of the smallest subsample.
4. Adjusts the weights when pooling cycles.
5. Uses strata and PSUs.
6. Uses domain estimation for subpopulations.
7. Reports the unweighted sample size (a STROBE reporting item, not in the NCHS text).
8. Handles missing data (a STROBE reporting item, not in the NCHS text).
9. Handles laboratory changes or calibration between cycles (bridging equations).
10. States the question addressed and defines the denominator, taking skip patterns into account (from the NCHS guidelines).
11. Considers statistical reliability (the NCHS presentation standards) when drawing inference in small groups (from the NCHS guidelines).
12. Does not estimate numbers of people by summing the final weights (from the NCHS guidelines; usually "not applicable").

Items 2 to 6 and 9 correspond to groups in the NCHS guidelines (weights, design variables, laboratory bridging and calibration); items 10 to 12 come from the NCHS guidelines of 26 January 2026 ([Guidelines for High Quality Analyses of NHANES Data](https://wwwn.cdc.gov/nchs/NHANES/QualityAnalysesGuidelines.aspx)). The NCHS item "read the NHANES documentation" cannot be rated from a paper and is not included.

**E3b is not an error rate.** In the abstract and title of the paper, E3b must not be presented as a proportion of wrong papers; the summary sentence says only that "most papers adjust for mediators in their association estimate, and this shrinks the OR by this much".

**Multiplicity.** No hypothesis tests. One primary estimate (E3a(A), primary denominator, unflagged group, classification "firm positive"). E3a(B) has the pre-declared reading rule above. Every other analysis is labelled "sensitivity, descriptive only" and is not used to change the conclusion.

Estimation: stratified weighted proportion, stratified variance with finite population correction, CI on the logit scale; if a stratum has 0% or 100%, use a stratified bootstrap by paper (2000 draws, percentile interval). Also report the unweighted proportion with a Wilson CI, and results by stratum (descriptive only).

Prior expectations (not decision thresholds): E3a(A) about 5 to 20%; E3a(B) higher than E3a(A), about 15 to 35%; E3a-label(A) about 20 to 40%; E3b above 50% with the OR shrinking toward 1; E1 denominator (i) about 50 to 70%; E2 about 20 to 40%; E4 about 3 to 6 threshold papers, most "not localisable" or inconclusive.

## 8. Why the criteria do not come out as expected by construction

| Criterion | What result would make it come out the other way | Control showing it can come out the other way |
|---|---|---|
| E3a | A paper that has all the G variables (DAG A), or is missing variables but adding them does not change the OR by more than 10% with certainty, is negative | The comparison only adds variables and drops none, so over-adjustment cannot push E3a up; the common sample removes the difference due to changing the sample; "firm positive" needs the CI outside the threshold; the negative control (random variable) is at most 0.05 and the positive control (drop then add age) shows the comparison can be positive when there is real confounding. The DAG A set of D1 contains only six social and behavioural variables that many papers already have, so a low E3a is a realistic outcome |
| E1 | A paper that describes enough and does it correctly reproduces | Three denominators; first combination locked in advance. Design information (pilot, data already seen): 3/4 papers passed |
| E2 | A paper that did the design correctly has a width ratio around 1 | Computed only on "near pass" papers, adjusted for N; the negative control is the distribution of the ratio across interpretation combinations plus rounding error; the 0.80 threshold is kept when fewer than 5% of values in this distribution fall below 0.80, otherwise it is lowered to the 5th percentile and recorded as a deviation |
| E4 | A true threshold of the same size is located; with no threshold there is no false localisation | Positive control at 50% slope change, conditions on coverage and false localisation |

The E3 of version 4 (outside the reference CI) and the definition of E3a through the step from step 3 to step 4 of the decomposition chain in version 5 (which merged adding confounders with dropping mediators) are dropped.

## 9. Conclusion rules and stopping (kill) criteria

### 9.1 Conclusion rule for each outcome (wording pre-declared)
Every result is reported because it is pre-registered.

| Outcome | "There is a problem" | "No evidence seen" | "Inconclusive" |
|---|---|---|---|
| E3a(A) | Lower bound of the 95% CI > 0.10 | Point estimate ≤ 0.10 and upper bound ≤ 0.25 | Everything else |
| E1 (denominator i) | Upper bound of the 95% CI < 0.80 | Lower bound ≥ 0.80 | Everything else |
| E2 | Lower bound of the 95% CI > 0.10 | Point estimate ≤ 0.10 and upper bound ≤ 0.25 | Everything else |

Pre-declared conclusion sentences for E3a:
- "There is a problem": "In this frame, at least about [lower bound] of papers omit confounders that are definite enough that adding them changes the OR by more than 10%."
- "No evidence seen": "No evidence was seen that the papers in the frame omit confounders definite enough to change the effect size; precision is limited (upper bound [x])." It must not be written as "the literature does not lack confounders".
- "Inconclusive": "The data are not sufficient to conclude on missing confounders."

E1 and E2 use similar sentences with the appropriate subject. The conclusion on each outcome stands alone; there is no rule of the type "if any outcome has a problem, the literature has a problem".

### 9.2 Stopping criteria
- **KD1, feasibility:** assessed after the specification table is locked and the indices are reconstructed for every paper, before any model is run. If fewer than 30 papers are testable: stop the full paper and write only a research letter on the TyG cluster.
- **KD2, extraction verification:** stop if any of the following fails after one prompt revision on the tuning set: error injection (Section 3.2, step 5); the test set, with the second test using the reserve test set (step 6); the author audit of the covariate list field, after one re-extraction (step 10).
- **KD3, implementation:** the variance implementation check (step 9) fails: stop until fixed; M1 to M4 are not run before that.
- **KD4, negative result for the main question:** if E3a(A) and E3a(B) both fall in "no evidence seen", the paper is a negative-result report on confounder specification (the main frame of the paper does not change into a story about reproduction or design after the numbers are seen); E1 and E2 are reported as secondary outcomes with their own conclusion sentences.
- **KD5, competing work:** before running M1 and before submission, scan Europe PMC and Semantic Scholar on three lines: (a) reproduction and survey-design checks of NHANES; (b) **assessment of the covariate sets of published literature against a DAG or for over-adjustment, with a reanalysis of effect sizes**; (c) computational reproduction on open data. If work already exists that does (b) on at least 10 NHANES index × kidney papers: drop the "confounder specification" framing, demote E3a to a secondary outcome and cite that work. If work already exists that does both reproduction and design checks on at least 10 NHANES index × kidney papers: E1 and E2 become descriptive only. Scan results are recorded in column (g) of the gap matrix.
- **KD6, time:** on working day 15 after registration, if the number of papers that have completed M1 is below 25 (based only on this count, with the time recorded): drop stratum F3. M1 results of F3 papers already completed are still reported, descriptively. The target estimand becomes "the eligible papers in the remaining strata", the weights are recomputed, and the narrowing is recorded as a pre-specified deviation.
- **KD7, enough papers for the main question:** if the primary denominator of E3a is below 20 papers, E3a is reported only as a count x/y (no proportion, no CI) and the paper is written as a descriptive research letter.

## 10. Transparency, deviations, contact with authors
- Everything open: source code, the frame file, the two screening and extraction passes and the reconciled version, the automatic check logs, error injection and author audit, the locked specification table and covariate label table, the log of sha256 of CDC files; archived with a DOI.
- Every deviation is recorded with date, reason and effect. Analyses conceived after seeing the data are labelled exploratory.
- The tone assesses methods and does not assign personal blame. By default, papers are named and the corresponding author is emailed at least 4 weeks before submission, with the results for their paper. If the author sends source code or clarifies the specification, rerun and record it as a deviation, reporting both numbers.
- The pilot part is reported as a pilot, separate from the main results.
- The search log for the gap matrix (queries, databases, dates) is published together with the matrix.

## 11. Known limitations
The frame comes only from Europe PMC and is filtered by title; results are not extrapolated beyond the frame. The DAGs are written by index family, not for each index; many nodes are "contested" and are handled with two DAGs, not with a single answer. Every paper is cross-sectional; reverse causation cannot be resolved. "Firm positive" is a conservative criterion, so E3a(A) may be lower than the true proportion. With the expected sample size, a negative conclusion for E3a only reaches "no evidence seen". There is no second human rater; extraction accuracy rests on error injection and the author's audit of three fields. The 0.80 threshold of E2 is not sensitive to small design factors. E4 rests on few papers.

## 12. Appendix A: DAG

The DAG figures and dagitty code for the six families D1, D2, D3, D4a, D4b, D5, under DAG A and DAG B, generated by source code from table 5.2, are attached to the OSF version (with the same hashes). Table 5.2 is the only source; if a figure and the table differ, the table is correct.

Integrity of the DAG folder (Edit B): the DAG files are listed with their sha256 in `dag/SHA256SUMS_i139_dag.txt` (25 files: 12 dagitty files, 6 png, 6 svg, and the manifest). The sha256 of that list file is `6a11943db20b5de0fdb1ea73bcc41257258e3a23756764fe1bcd125ab091aacb`. The manifest `dag_manifest_i139.json`, which records the code check that for all 12 DAGs the pre-declared adjustment set satisfies the backdoor criterion and is minimal, has sha256 `7dc7cfc6b1f8f2b6cde66782bb12c46f6882f1465a6755da803448a4e9159abc`. The DAGs are drawn at the level of index families and only edges that decide the adjustment set are drawn, so the backdoor check shows that the code matches the table, not that the table is correct.

## 13. Response to the version 5 review (round 2)

| # | Reviewer point | Decision | Change |
|---|---|---|---|
| 1 (major) | E3a has no denominator; the comparison from step 3 to step 4 merges adding confounders with dropping mediators | Partly accepted | A separate comparison "paper's set + missing variables" versus "paper's set", on the common sample; two-order decomposition (M3). The primary denominator is "testable" papers without the flag, not "near pass" papers as the reviewer proposed: E3a asks what the paper's own specification is missing, which does not depend on reproducing the exact number; "near pass" is secondary denominator (i), with lower and upper bounds for papers that are not testable (Section 7) |
| 2 (major) | False positives from changing the sample and noise; the OR is non-collapsible | Accepted | Common sample, with the sample-change difference reported separately; three classification levels, main E3a uses "firm positive"; PR and marginal standardisation sensitivity analyses; added negative control (random variable) and positive control (age) for E3a (M3, Section 8) |
| 3 (major) | No variable equivalence rule | Accepted | Table 5.4, labels assigned by source code |
| 4 (major) | DAG lacks nodes, contested ranks are forced; F1 to F4 are inconsistent on diabetes; F5 left for later | Accepted | Full node table by six DAG families with ranks G/T/K/P and a one-line reason (5.2); DAG A and DAG B reported in parallel with a reading rule; diagnosed diabetes in D4a moved to "contested, in DAG B", with a justification different from D1 (D1 contains glucose); D5 has its own column, family assigned by source code (5.3); reverse causation recorded in the DAG. Figures and dagitty code generated from the table in Appendix A |
| 5 (major) | Extraction verification is not enough for the fields that decide E3a | Mostly accepted | 60 errors, half "plausible but wrong", every error keeps a real excerpt; condition on the CI lower bound at least 0.85; two model families for the E3a fields, and if not available, agreement between passes is not used as evidence; κ and Gwet AC1 by field; author audit made mandatory. Partly accepted on the audit scope: only three fields, for papers with a "missing" label and 10 random papers, not for every paper entering M3, because a paper missing no variable counts as negative without a model, and for "near pass" papers reproduction is already objective evidence of extraction |
| 6 (major) | KD4 almost never triggers; "any one" type of conclusion | Accepted | Conclusion rule for each outcome based on CIs, wording pre-declared (9.1); new KD4 for a negative result on the main question; KD7 when the E3a denominator is below 20 |
| 7 (major) | Scope does not match the TyG and confounding question | Accepted | The target estimand is the whole frame; TyG reported as x/y; DAG for every family at the same level of detail (5.2); E3a is the only primary outcome, E1/E2 become secondary outcomes |
| 8 | Test set reused after prompt revision | Accepted | Reserve test set of 5 papers locked in advance (3.2 step 6, KD2) |
| 9 | Gap matrix lacks a confounding/DAG column | Partly accepted | KD5 gains line (b) and the gap matrix gains column (g). Done 14 Nov 2026; KD5 row (b) not triggered (see Related work) |
| 10 | The Kish example does not match | Accepted | Recomputed with F1 = 7.5 papers: factor 1.121, effective size 42.35, half-width 0.1505 (0.1259 with the finite population correction) (2.3) |
| 11 | 1.563 should be 1.562 | Accepted | Section 7 |
| 12 | KD6 may be influenced by results already seen | Accepted | Based only on the number of papers that completed M1 on day 15; F3 results already available are still reported |
| 13 | The frame query is not in the registration | Accepted | Copied verbatim into 2.1 |
| 14 | Multiplicity | Accepted | One primary estimate; everything else "sensitivity, descriptive only" (Section 7) |
| 15 | Restriction to people with diabetes; reverse causation | Accepted | Label "conditioning on a descendant of the exposure" (2.2); reverse-causation paragraph in 5.2 |
| 16 | A sentence forbidding E3b as an error rate | Accepted | Section 7 |

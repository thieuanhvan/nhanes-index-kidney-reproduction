# Dietary index reconstruction notes, I139 (outcome-blind stage)

Code: tools/diet_indices_i139.py (DIET dict keyed by rec_id; reading 1 / reading 2 per paper). Output: extraction/diet_recon_i139.tsv.
Built 2026-11-15 on the outcome-blind table (nhanes_indices.load_blind: no creatinine, urine albumin, eGFR or CKD column is read). No outcome model fitted.
Numbers are unweighted, among participants aged >= 18 with a non-missing index, truncated to 4 significant digits.

Common rules
- Recall used only when its status is 1 (reliable, met minimum criteria). "Mean of 2 days" = mean of day 1 and day 2 when day 2 is valid, else day 1. 1999-2002 have one recall only.
- Tertile cut points: 1/3 and 2/3 sample quantiles (linear interpolation) in the stated reference sample; value <= q1 is tertile 1, <= q2 tertile 2, else tertile 3. Ties at a cut point (zero alcohol) go to the lower tertile.
- Physical activity (2007-2018): GPAQ MET-min/week, MET 8 vigorous work and leisure, 4 moderate work, walk/bicycle and moderate leisure; all five domains required for the total.
- Physical activity (1999-2006, update 2026-11-15): individual activity files PAQIAF, PAQIAF_B, PAQIAF_C, PAQIAF_D (added to nhanes_extra_i139.STEMS as "PAQIAF"; read directly because there is one row per activity). Leisure MET-min/week = sum over activities of MET x times in past 30 days x minutes x 7/30; 0 when both PAD200 and PAD320 are answered no/unable; missing when any activity row lacks times or minutes. MET is either 8/4 by PADLEVEL (vigorous/moderate, the GPAQ weights) or the activity's own PADMETS, as each paper defines. Transport: PAD020 walk/bicycle, PAQ050Q per day/week/month (PAQ050U) converted to times per week, x PAD080 minutes x 4 MET. Work minutes were not collected in 1999-2006, so the 1999-2006 total is leisure + transport and runs lower than the 2007-2018 GPAQ total (adult medians about 210-420 vs 960-1320 MET-min/week); pooled tertiles therefore put more 1999-2006 people in the lowest tertile. Adults with a PA value: about 5,400-6,000 per cycle in every cycle.
- Vitamin D (DR1TVD) exists only from 2007-2008. 1999-2000 has no individual carotenoids (only DRXTCARO, carotene in RE) and vitamin E only as DRXTVE (alpha-tocopherol equivalents); 2001+ uses alpha-tocopherol (ATOC).
- Cotinine 2005-2006 from COT_D (coordinator update 2026-11-15).
- No energy adjustment: none of the 7 papers states one for the index.
- DII constants: Shivappa N et al., Public Health Nutr 2014;17(8):1689-1696, doi [[DOI withheld]](https://doi.org/[DOI withheld]), copied from the open-source R package dietaryindex, [R/DII.R](https://github.com/jamesjiadazhan/dietaryindex/blob/main/R/DII.R) (last commit on the file f478683d1e7130b731c53ad5b078f39481051161, 2024-09-29; file sha256 c6392632...dc2a; fetched 2026-11-15). NHANES column mapping as in the package (R/DII_NHANES_FPED.R, nutrient part only, no FPED): n-3 = P183+P184+P205+P225+P226, n-6 = P182+P204, vitamin A = VARA (RAE), vitamin E = ATOC, folic acid = FA, caffeine mg/1000. Score per parameter = effect x (2 Phi(z) - 1).

## 71eedfeafb (OBS, 1999-2018)

Definition used. Quote (Methods, Assessment of OBS, para 1): "Each of the 20 OBS components was categorized into sex-specific tertiles." Also: "[excerpt withheld]" and "Participants were required to have valid data for at least 19 of the 20 OBS components". 16 nutrients (fiber, carotenoids, riboflavin, niacin, B6, B12, C, E, total folate, calcium, magnesium, zinc, copper, selenium; total fat and iron pro-oxidant) plus physical activity (MET), alcohol (from the recalls), serum cotinine, BMI (pro-oxidant). Scores 0-2, sum of available scores when >= 19 present.

Decisions
- Tertile reference sample. reading_1: the paper's population without outcome criteria (age >= 20, BMI >= 25, not pregnant). reading_2: all participants of the pooled cycles.
- "Carotenoids". reading_1: sum of alpha-carotene, beta-carotene, beta-cryptoxanthin, lycopene, lutein+zeaxanthin. reading_2: beta-carotene only.
- Alcohol. reading_1: sex-specific tertiles of recall alcohol g/d, reversed (literal "each of the 20"; zeros fall in tertile 1). reading_2: 2/1/0 non / non-heavy / heavy drinker.
- Physical activity: tertiles of total MET-min/week (work, transport, leisure).

- Physical activity 1999-2006 (update): leisure from PAQIAF with MET 8 vigorous / 4 moderate by PADLEVEL (the paper's "[excerpt withheld]"), plus walk/bicycle transport at 4.0.

Components unavailable: 1999-2000 carotenoids (no individual carotenoids in DRXTOT). The 19-of-20 rule is applied as written, so 1999-2000 participants are scored on the other 19 components when all 19 are present: 4,310 adults (>= 18) scored in 1999-2000, 2,461 of them in the outcome-free population.

Status: computed. n_index 51,920 (reading_1 mean 20.77, SD 7.354; reading_2 mean 21.19). In the outcome-free population (age >= 20, BMI >= 25, not pregnant) n 32,648, mean 20.34; the paper reports 19.84 (weighted, n 28,377).

## 9e2e6087fa (DII, 2005-2018)

Definition used. Quote (Materials and methods, Definition of dietary inflammation index, para 1): "[excerpt withheld]". The list names 27 parameters (checked in the XML too); vitamin B6 is the one of the usual 28 NHANES parameters that is not named. Day: "consumed by each participant during the day (midnight to midnight)".

Decisions
- Components. reading_1: the 27 named parameters. reading_2: the 28 (adds vitamin B6).
- Days. reading_1: day-1 recall ("the day"). reading_2: mean of the two recalls (intakes averaged, then scored). In the code reading_2 changes both points together.
- Prereg section 4 applied strictly (coordinator decision 2026-11-15): the paper says 28 components but names 27, so the index is not reconstructable. The 27-component version is kept only as a sensitivity row (status sensitivity_only), with the B6 + 2-day version as its reading_2.

Components unavailable: 2005-2006 vitamin D (parameter skipped for that cycle, as dietaryindex does).

Status: not_reconstructable (paper says 28 components but names 27; vitamin B6 not named). Sensitivity: n 37,483, mean 1.254, SD 1.804; its reading_2 mean 1.066, SD 1.916.

## 988df86a79 (DII, 2011-2018)

Definition used. Quote (Methods, Definition of the DII and CKD, para 1): "[excerpt withheld]" Components: energy, carbohydrate, protein, total fat, fiber, cholesterol, SFA, MUFA, PUFA, n-3, n-6, vitamins A, B1, B2, niacin, B6, B12, C, D, E, folic acid, alcohol, beta-carotene, caffeine, iron, magnesium, zinc, selenium (28). Formula quoted: "Zscore1 = Zscore → (converted to percentile)× 2 − 1"; constants "provided in the study by Shivappa N et al."

Decisions: reading_1 day-1 recall (stated); reading_2 mean of two recalls (sensitivity only).

Components unavailable: none.

Status: computed. n_index 20,681, mean 1.107, SD 2.001 (median 1.328). The paper reports mean 1.50 (SD 1.91), median 1.74, n 19,317. The gap of about 0.4 is not explained by the definition as written. Checked against the dietaryindex mapping (update 2026-11-15), adults, day 1: folic acid DR1TFA (package, kept) gives 1.107; total folate DR1TFOLA would give 0.942; DFE (DR1TFDFE) 0.898; caffeine left in mg instead of g 0.962. Vitamin A: 2011-2018 files carry only RAE (DR1TVARA, used by the package) and retinol, no IU or RE column. Beta-carotene (DR1TBCAR, ug) and alcohol (DR1TALCO, g) match the package and Shivappa units. No column or unit is wrong by the paper's or the package's definition, and every alternative moves the mean further from 1.50, so nothing was changed. Largest positive contributions: beta-carotene +0.35, vitamin D +0.22, vitamin C +0.20, vitamin A +0.20, alcohol +0.18, fiber +0.18. Left for M1 (sample and weighting may explain it).

## fe82b9da6a (OBS, 2011-2014)

Definition. Quote (Methods, Main variables): "[excerpt withheld]" Lifestyle: physical activity (MET), alcohol (nondrinker / non-heavy / heavy, men 30 g, women 15 g), BMI, cotinine. "All factors were grouped based on tertiles." The 16 dietary components are never named (text, tables, XML; the only supplement is a figure); the method is cited to ref 20 (Ke et al., Front Physiol 2023).

Decisions
- reading_1: not reconstructable (prereg section 4: the paper does not name all components).
- reading_2 (sensitivity only, not the primary status; status sensitivity_only row in the tsv): common 16-nutrient OBS list as named by aca6f864fb, day-1 recall, tertiles not sex-specific among diabetes_any adults >= 20, alcohol scored 2/1/0, all 20 required.

Components unavailable: none for reading_2.

Status: not_reconstructable (16 dietary components unnamed). Sensitivity reading_2 n 9,787, mean 22.16, SD 7.206.

## aca6f864fb (OBS, 2005-2018)

Definition used. Quote (2.3 OBS): "[excerpt withheld]" Alcohol 2/1/0 (non, moderate women 0-15 g/d men 0-30 g/d, heavy); BMI obesity 0, overweight 1, normal 2; physical activity = MET x frequency x duration per week; pro-oxidants total fat, iron, cotinine, alcohol, obesity. Exclusion "items comprising the OBS were < 16" read as: score kept when at least 16 of the 20 components are present.

Decisions
- Days. reading_1: mean of the recalls ("24-h dietary recall interviews", plural). reading_2: day 1.
- Tertile sample. reading_1: adults >= 18 (the paper's age limit). reading_2: all participants.
- "Carotene". reading_1: alpha + beta carotene. reading_2: beta-carotene only.
- "Total folic acid" read as total folate (DR1TFOLA) in both readings.

- Physical activity 2005-2006 (update): PAQIAF leisure with each activity's own PADMETS ("metabolic equivalent of each physical activity × frequency of physical activity × duration") plus walk/bicycle transport at 4 MET; 2007-2018 GPAQ total.

Components unavailable: none (2005-2006 physical activity restored from PAQIAF).

Status: computed. n_index 37,483, mean 20.37, SD 7.240 (reading_2 mean 21.49) (paper quartile limits [3,15], [15,20], [20,26], [26,37]).

## 9835a67923 (DOBS, 2001-2010)

Definition used. Quote (Materials and methods, Evaluation of DOBS): "[excerpt withheld]" and "The individual scores for these 16 components were summed to derive the DOBS". Components from Table 1 rows (fiber, carotene, B2, niacin, B6, total folate, B12, C, E, calcium, magnesium, zinc, copper, selenium, total fat, iron); pro-oxidants fat and iron. The scoring rule is in Supplementary Tables 1-2 (PMC supplement Table_1.DOCX, not downloaded); the main text gives none.

Decisions
- Scoring: tertiles 0-2, reversed for fat and iron (the paper's mean 16.81 on a 0-32 range fits this).
- Sex-specific. reading_1: not sex-specific (not stated). reading_2: sex-specific.
- Tertile sample: outcome-free population (diabetes_any, not pregnant, energy 500-5000 kcal/day) in both readings.
- "Carotene". reading_1: alpha + beta. reading_2: beta only.
- All 16 required (paper excludes missing dietary intake).

Components unavailable: none (2001-2010 all have the 16).

Status: computed. n_index 27,054, mean 16.93, SD 7.096 (paper: 16.81 +/- 6.61, diabetes sample). The supplement would settle the sex-specific point. Download attempt 2026-11-15 (coordinator-approved, once): the PMC URL (pmc.ncbi.nlm.nih.gov/articles/instance/11965115/bin/Table_1.DOCX) returned a proof-of-work bot challenge page and the Europe PMC URL (europepmc.org/articles/PMC11965115/bin/Table_1.DOCX) returned HTTP 403. Both are blocks; stopped, nothing saved, reading_1 unchanged.

## b966d3cfdc (OBS, 2007-2018)

Definition used. Quote (Materials and methods, Oxidative balance score, p. 7): "Antioxidants received scores of 1, 2, and 3 from lowest to highest tertiles"; "Table 5 provides the categorization and scoring criteria for each component included in the OBS." Table 5 (p. 8) prints sex-specific cut points for all 16 non-alcohol components; these fixed cut points are used directly (no tertiles recomputed). Antioxidants: vitamin C, vitamin E, alpha-carotene, beta-carotene, beta-cryptoxanthin, lycopene, lutein+zeaxanthin, selenium, zinc, physical activity. Pro-oxidants: total fat, PUFA, SFA, iron, BMI, cotinine, alcohol (3 none, 2 non-heavy, 1 heavy). Interval rule from the table: [a, b) is the middle score.

Decisions
- Days. reading_1: day-1 recall ("a 24-hour dietary recall interview"). reading_2: mean of the recalls.
- Physical activity. reading_1: leisure-time MET-min/week ("calculated based on leisure-time activity data"). reading_2: all GPAQ domains.
- All 17 required (paper excludes incomplete OBS data).

Components unavailable: none.

Status: computed. n_index 30,415, mean 34.73, SD 4.384 (diabetes_any adults: mean 34.40).

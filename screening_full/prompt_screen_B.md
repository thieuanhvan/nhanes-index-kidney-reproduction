# Screening prompt B: I139 eligibility check from title and abstract

## Task in brief

Read `screening_full/packet_screen.tsv`. It is a tab-separated UTF-8 file with one row per paper and the columns `rec_id`, `title`, `abstract`, `year`. For every row, record whether the paper passes the abstract-level eligibility check of study I139, why, two yes/no labels, and a supporting excerpt. Work from the packet text alone: no internet lookup, no recollection of the paper, and no access to any other screening result in the folder.

## What to deliver

A tab-separated UTF-8 file, saved where the operator asks (for example `screening_full/screen_B.tsv`), with exactly this first line:

```
rec_id	decision	reason_code	flag_outcome_component	flag_exposure_descendant	quote
```

followed by exactly one line for each packet row, in the same order as the packet, and nothing else.

- `decision` takes one of three lower-case values: `include`, `exclude`, `unclear`.
- `reason_code` must belong to the decision: `INC_ALL_MET` for `include`; a code from table 1 for `exclude`; a code from table 2 for `unclear`.
- The two flag columns hold `0` or `1`.
- `quote` holds a single passage taken word for word from the title or the abstract (see step 5). It must not contain a tab or a line break and is not wrapped in quotation marks.

## Step 1. Check the grounds for exclusion

Go through table 1 from the top. Stop at the first row that the title or abstract positively shows to be true, and use its code with `decision` = `exclude`. Silence in the abstract is never a ground for exclusion.

Table 1. Exclusion codes, in the order they are checked

| Order | Code | The title or abstract shows that |
|---|---|---|
| 1 | `EX_PUB_TYPE` | the record is a comment, letter, editorial, erratum, a review of any kind including systematic reviews and meta-analyses, or a pooled analysis across studies |
| 2 | `EX_RETRACTED` | the paper has been retracted |
| 3 | `EX_NOT_NHANES` | NHANES is not the main data source, or the NHANES data used are only NHANES III (1988 to 1994) or only fall outside the 1999 to 2023 cycles |
| 4 | `EX_KIDNEY_EXPOSURE` | the exposure is UACR, eGFR or another measure of kidney function or damage, not the outcome |
| 5 | `EX_NOT_INDEX` | the main exposure is a single quantity or behaviour (one nutrient, one pollutant level, hours of sleep and similar) and not an index, ratio or score |
| 6 | `EX_RESTRICTED_DATA` | the index can only be computed with restricted NCHS data, such as geocoded or linked area-level information or restricted-use files |
| 7 | `EX_CKD_MORTALITY_PROGRESSION` | the outcome is death or progression of kidney disease among people who already have CKD |
| 8 | `EX_OUTCOME_NOT_KIDNEY` | no outcome is prevalent kidney disease (CKD, albuminuria, low eGFR, or DKD) |
| 9 | `EX_NOT_CROSS_SECTIONAL` | the analysis is only longitudinal or time-to-event, and no cross-sectional part can be separated |
| 10 | `EX_NO_ADJUSTED_ESTIMATE` | the paper gives no adjusted association estimate, for example it reports only prediction performance, only prevalence figures or only crude comparisons |

A paper with several data sources is not excluded under code 3 if its NHANES part can be analysed on its own. A title containing "albumin-to-creatinine ratio" as the outcome does not make the exposure an index.

## Step 2. Check the five eligibility conditions

If no row of table 1 applies, check that the title and abstract support all of the following. These are the study's eligibility conditions with the reading rules used at abstract level.

(a) The paper reports an adjusted association estimate (an OR, PR or beta with a confidence interval or standard error). At abstract level, a reported adjusted OR, PR or beta, or a statement that multivariable or adjusted regression (logistic, linear or similar) was fitted, is enough.

(b) The outcome is kidney disease present at the time of the survey: CKD, albuminuria, low eGFR, or CKD in people with diabetes (DKD, diabetic nephropathy). One such outcome among several is enough.

(c) The main exposure is an index, ratio or score computable from public NHANES variables, such as a diet quality or dietary inflammation index, a blood-count inflammation ratio, a lipid ratio or a triglyceride-glucose index.

(d) The data come mainly from NHANES, continuous cycles 1999 to 2023. "NHANES" with no years stated is accepted.

(e) The analysis is cross-sectional or has a separable cross-sectional part. Exposure and kidney outcome measured at the same NHANES examination count as cross-sectional, whether or not the abstract uses the word.

If all five are supported, set `decision` = `include` and `reason_code` = `INC_ALL_MET`.

## Step 3. Otherwise mark the record as unclear

If nothing in table 1 applies but one or more of the five conditions cannot be judged from the text, set `decision` = `unclear`. Pick the code from table 2 for the first condition, taking them in the order of the table, that cannot be judged.

Table 2. Unclear codes, in the order they are checked

| Order | Code | What cannot be judged |
|---|---|---|
| 1 | `UN_DATA_SOURCE` | whether NHANES is the main data source |
| 2 | `UN_EXPOSURE` | what the main exposure is, or whether it is an index, ratio or score from public NHANES variables |
| 3 | `UN_OUTCOME` | whether the outcome is prevalent CKD, albuminuria, low eGFR or DKD; this also covers a kidney outcome that is only eGFR or UACR as a continuous variable |
| 4 | `UN_DESIGN` | whether a cross-sectional analysis exists |
| 5 | `UN_ESTIMATE` | whether any regression or adjusted estimate is reported, because the abstract mentions none |

When the abstract field is empty and the title does not show a ground in table 1, use `UN_NO_ABSTRACT` instead of the codes in table 2.

## Step 4. Set the two labels on every record

Fill both columns for all records, including excluded ones. Look at the main index, or at any of them if the paper has more than one main index. Take its components from the abstract; if the abstract gives only the name of a well-known index, use that index's standard published formula. A record with no index gets 0 in both columns.

`flag_exposure_descendant` records the descriptive label "conditioning on a descendant of the exposure". It is 1 when both hold: the analysis is limited to people with diabetes (for instance DKD among diabetic adults), and the index includes glucose, insulin or HbA1c (TyG and its BMI, waist or waist-to-height versions, HOMA-IR, METS-IR and the like). In every other case it is 0.

`flag_outcome_component` is 1 when the index includes something that the outcome itself determines. The list is: serum albumin, haemoglobin, creatinine, cystatin C, blood urea nitrogen, urinary albumin, potassium, phosphorus, and uric acid only if the outcome involves eGFR. A CKD outcome defined through eGFR, alone or together with albuminuria, involves eGFR; an albuminuria-only outcome does not. Indices that set this label include the uric acid to HDL cholesterol ratio (when the outcome involves eGFR), the BUN to creatinine ratio, the advanced lung cancer inflammation index, the geriatric nutritional risk index, the prognostic nutritional index, the neutrophil percentage to albumin ratio, the C-reactive protein to albumin ratio, and the HALP score (haemoglobin, albumin, lymphocytes, platelets). Otherwise it is 0.

Neither label changes the decision.

## Step 5. Copy the supporting excerpt

Put in `quote` the passage of the title or abstract that best supports the decision and the code, copied character for character (same case, spelling and punctuation). It must be one unbroken stretch, no longer than 40 words, with no ellipsis, no merging of separate passages and no rewording. For `UN_NO_ABSTRACT`, copy the whole title.

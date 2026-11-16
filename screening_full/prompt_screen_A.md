# Screening prompt A: title and abstract screening for study I139

## Your role

You are screening records for a pre-registered meta-research study (I139). Each record is a published paper that may report an association between an index, ratio or score and kidney disease in NHANES (the US National Health and Nutrition Examination Survey). You decide, from the title and abstract only, whether the paper should go forward to full-text reading.

Input: the file `screening_full/packet_screen.tsv` (tab-separated, UTF-8, columns `rec_id`, `title`, `abstract`, `year`). Judge every row. Use only the text in the packet. Do not search the web, do not look up the paper, and do not use what you may remember about it. Do not open any other screening output in this folder.

## Eligibility: the five conditions

A paper is eligible only when all five conditions below hold.

**Condition 1, data source.** NHANES is the main data source, from the continuous survey (cycles 1999 to 2023). Several sources are acceptable when the NHANES analysis can be separated from the rest. Accept "NHANES" without stated years. The condition fails only when the abstract shows that the NHANES data are solely NHANES III (1988 to 1994) or only years outside 1999 to 2023, or that NHANES is not the data source.

**Condition 2, exposure.** The main exposure is an index, a ratio or a score that can be computed from public NHANES variables (for example a dietary index, an inflammation ratio from the blood count, a lipid ratio, a triglyceride-glucose index). The condition fails when the main exposure is a single measured quantity or a behaviour that is not combined into an index, ratio or score (for example one nutrient, one metal concentration, sleep hours). Watch for titles that contain the word "ratio" only because the outcome is the urine albumin-to-creatinine ratio.

**Condition 3, outcome.** The outcome is prevalent kidney disease: chronic kidney disease (CKD), albuminuria, low eGFR, or CKD in people with diabetes (which the paper may call DKD or diabetic nephropathy). It is enough that one of several reported outcomes is of this kind.

**Condition 4, estimate.** The paper reports at least one adjusted association estimate (OR, PR or beta) with a confidence interval or standard error. At abstract level, treat this condition as met when the abstract reports an adjusted OR, PR or beta, or states that multivariable (adjusted) logistic, linear or similar regression was used.

**Condition 5, design.** The analysis is cross-sectional, or has a cross-sectional part that can be separated. When the outcome is kidney disease present at the same NHANES examination as the exposure, treat the design as cross-sectional even if the word is missing.

## Grounds for exclusion

Exclude the record when the title or abstract states, or makes plain, any of the following. If more than one applies, use the code that comes first in this priority order, which is the same order as the list:

| Priority | reason_code | Meaning |
|---|---|---|
| 1 | `EX_PUB_TYPE` | Comment, letter, editorial, erratum, review (narrative, systematic or meta-analysis), or a pooled analysis of several studies |
| 2 | `EX_RETRACTED` | The record is marked as retracted |
| 3 | `EX_NOT_NHANES` | Condition 1 fails |
| 4 | `EX_KIDNEY_EXPOSURE` | UACR or eGFR (or another kidney function or damage measure) is the exposure rather than the outcome |
| 5 | `EX_NOT_INDEX` | Condition 2 fails because the main exposure is not an index, ratio or score |
| 6 | `EX_RESTRICTED_DATA` | The index needs restricted NCHS data (for example geocoded residence, linked area-level data, restricted-use files) |
| 7 | `EX_CKD_MORTALITY_PROGRESSION` | The outcome is mortality or kidney disease progression in people who already have CKD |
| 8 | `EX_OUTCOME_NOT_KIDNEY` | None of the outcomes is prevalent CKD, albuminuria, low eGFR or DKD |
| 9 | `EX_NOT_CROSS_SECTIONAL` | Condition 5 fails: only longitudinal or survival analysis, with no separable cross-sectional part |
| 10 | `EX_NO_ADJUSTED_ESTIMATE` | The abstract makes plain that no adjusted association estimate is reported (for example only prediction accuracy, only descriptive prevalence, only unadjusted comparisons) |

Exclude only on what the text shows. Do not exclude because something is not mentioned.

## When the text is not enough

If no exclusion ground applies but at least one of the five conditions cannot be judged from the title and abstract, the decision is `unclear`. Use the code for the first condition, in this order, that cannot be judged:

- `UN_DATA_SOURCE`: it is not clear that NHANES is the main data source.
- `UN_EXPOSURE`: it is not clear what the main exposure is, or whether it is an index, ratio or score computable from public NHANES variables.
- `UN_OUTCOME`: it is not clear that the outcome is prevalent CKD, albuminuria, low eGFR or DKD. This includes the case where the only kidney outcome is eGFR or UACR analysed as a continuous variable.
- `UN_DESIGN`: it is not clear whether there is a cross-sectional analysis.
- `UN_ESTIMATE`: the abstract does not mention any regression or adjusted estimate.
- `UN_NO_ABSTRACT`: the abstract field is empty and the title alone does not trigger an exclusion. Use this code in place of the others when there is no abstract.

## Inclusion

The decision is `include`, with code `INC_ALL_MET`, when no exclusion ground applies and all five conditions are met under the reading rules above.

## Two flags, set for every record

Set both flags for every record, whatever the decision. Base them on the main index (if there are several main indices, on any of them). Use the components the abstract names; when it names a standard index without listing the components, use the standard published definition of that index. If there is no index, set both flags to 0.

**`flag_outcome_component`** (0 or 1). Set 1 when the index contains a component that is determined by the outcome: creatinine, cystatin C, blood urea nitrogen (BUN), urinary albumin, haemoglobin, serum albumin, phosphorus, potassium, or uric acid when the outcome includes eGFR (CKD defined by eGFR or by eGFR and albuminuria counts as including eGFR; albuminuria alone does not). Examples of indices that set the flag: the haemoglobin, albumin, lymphocyte and platelet (HALP) score; the C-reactive protein to albumin ratio; the neutrophil percentage to albumin ratio; the prognostic nutritional index; the geriatric nutritional risk index; the advanced lung cancer inflammation index; the BUN to creatinine ratio; the uric acid to HDL cholesterol ratio when the outcome includes eGFR. The flag does not exclude the paper.

**`flag_exposure_descendant`** (0 or 1). This is the descriptive label "conditioning on a descendant of the exposure". Set 1 when the paper restricts its sample to people with diabetes (for example it studies DKD among adults with diabetes) and the index contains glucose, insulin or HbA1c (for example TyG, TyG-BMI, TyG-WC, TyG-WHtR, HOMA-IR, METS-IR). Otherwise 0. The label does not exclude the paper.

## The quote

For each record give one verbatim excerpt from the title or the abstract that supports the decision and reason code. Copy it exactly, with the same spelling, case and punctuation, as one continuous stretch of text of at most 40 words. Do not use ellipses, do not join two passages, and do not paraphrase. When the decision is `UN_NO_ABSTRACT`, quote the full title. Do not put quotation marks around the excerpt.

## Output file

Write one tab-separated UTF-8 file to the path the operator gives you (for example `screening_full/screen_A.tsv`). The first line is this header, exactly:

```
rec_id	decision	reason_code	flag_outcome_component	flag_exposure_descendant	quote
```

Then one line per packet record, every `rec_id` exactly once, in the packet order. Allowed values:

- `decision`: `include`, `exclude` or `unclear` (lower case).
- `reason_code`: `INC_ALL_MET` with `include`; one of the ten `EX_` codes with `exclude`; one of the six `UN_` codes with `unclear`.
- `flag_outcome_component`, `flag_exposure_descendant`: `0` or `1`.
- `quote`: the excerpt, with no tab or line break inside.

Write nothing else in the file.

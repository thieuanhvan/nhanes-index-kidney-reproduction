# NCHS/STROBE reporting checklist: brief (written 15/11/2026 before rating)

You rate how published NHANES papers report their analysis. Working folder:
<workdir>/idea-hunt-metaresearch-20261114

Read ONLY: `extraction/sample_47.tsv` (column rec_id gives the 47 papers) and `fulltext_txt/<rec_id>.txt` for each
(if a table or formula is missing, `fulltext/<rec_id>.xml`). Do not read other files; no web.

For each paper rate 12 items as `yes` / `no` / `not_applicable`, with the line number(s) of the evidence:
1. States the survey cycles.
2. States the name of the weight variable (e.g. WTMEC2YR, WTSAF2YR, WTDRD1); "weighted" without the name = no.
3. Uses the weight of the smallest subsample among the variables analysed (e.g. fasting weight when fasting glucose or
   TG or LDL is used; dietary weight for dietary recall indices). If the paper is unweighted, = no.
4. Adjusts the weights when pooling cycles (e.g. divides 2-year weights by the number of cycles, 4-year weights for
   1999-2002, 3.2 years for 2017-March 2020). Unweighted = no.
5. Uses strata and PSUs (SDMVSTRA, SDMVPSU, "complex survey design", svy commands with design variables).
6. Uses domain (subpopulation) estimation for subpopulations rather than subsetting the data; not_applicable if the
   analysis is the whole adult sample with only complete-case exclusions; if the paper subsets a subpopulation
   (e.g. diabetes only) and does not mention subpopulation/domain analysis = no.
7. Reports the unweighted sample size of the analytic sample.
8. States how missing data were handled.
9. Handles laboratory changes or calibration between cycles (e.g. creatinine calibration for 1999-2000 and 2005-2006,
   HbA1c or glucose bridging). not_applicable if all cycles are 2007 or later and no such variable crosses a method
   change mentioned in NHANES documentation; if the paper uses 1999-2000 or 2005-2006 serum creatinine and says nothing
   = no.
10. States the question addressed and defines the denominator, taking skip patterns into account (e.g. how
    questionnaire skips for non-drinkers or never-smokers were coded).
11. Considers statistical reliability in small groups (NCHS presentation standards, relative SE, suppression of
    unstable estimates). Usually no.
12. Does not estimate population numbers by summing final weights incorrectly; not_applicable if no population counts
    are reported; yes if counts are reported with correct method.

Output: write `extraction/nchs_checklist_i139.tsv` with columns `rec_id, item, rating, evidence_line` (12 rows per
paper, 564 rows). Then report only the counts of yes/no/not_applicable per item.

# Extraction prompt B (I139)

You are extracting structured facts from the full text of ONE published paper (NHANES, an index as exposure, kidney disease as outcome). The text is given to you as plain text. Your whole answer is a single JSON object and nothing else (no prose before or after). Everything you report must be traceable to the text. You report what the authors wrote, not what you think is correct.

## 1. Output format

Use exactly the schema below. Key names must match character for character.

A "field" is an object with three keys:
`{"value": <the value>, "quote": "<text copied from the paper>", "loc": "<where in the paper>"}`

- `quote` must be copied character for character from the full text. The text comes from a PDF and has line breaks inside words ("Albuminu-" at the end of one line, "ria" at the start of the next), ligatures, HTML entities such as `&lt;`, and odd spacing. Copy these exactly as they appear (a line break may be written as a single space), and do not repair them. Better still, choose a quote span that contains no hyphen-broken word and does not cross a page break (the "## PAGE n" headers and running page headings are not part of the sentence). Never paraphrase, never fix typos, never merge two separate places into one quote, never insert "..." inside a quote. If you need two pieces of text, choose the one piece that carries the value. Keep quotes short (one sentence, one table row, one footnote) but long enough to contain the value: the quote must itself show the value you put in `value` (a quote that only names the index does not support its formula).
- `loc` names the place: "Abstract", "Methods 2.3 para 2", "Table 2", "Table 2 footnote", "Results para 1", etc. When the value comes from a table, loc MUST contain the word Table and its number (e.g. "Table 3 footnote"). An automatic check uses this.
- If the paper does not report the item, use `{"value": null, "quote": "", "loc": ""}`. Do not guess and do not fill from background knowledge.
- Numbers go in `value` as JSON numbers, not strings. Lists go in `value` as JSON arrays of strings.

Top level keys:

- `rec_id`, `doi`: plain strings (rec_id is given to you; take the DOI from the text if present, else "").
- `index_name` (field), `index_formula` (field, formula as written), `index_components` (field, value = list of component names such as ["TG","FPG"]).
- `cycles` (field, value = list such as ["2011-2012","2013-2014"]; write each cycle as YYYY-YYYY).
- `population` (field: age limits and main inclusion/exclusion as written).
- `n_analytic` (field, integer: the number of people in the model that gives the main estimate. If that model uses only a subset (for example people above a turning point, or one sex), give the subset N and write the overall analytic N in `notes`. Otherwise give the final analytic sample).
- `outcomes_in_abstract` (field, value = list of outcomes named in the abstract).
- `outcome_main` (field: outcome of the main estimate, chosen by rule R1).
- `outcome_definition` (field), `egfr_equation` (field), `uacr_threshold` (field).
- `weight_variable` (field: name or description of the survey weight used), `cycle_pooling` (field: how weights were combined across cycles).
- `main_estimate` (object) with fields: `model_label` (e.g. "Model 3"), `effect_type` (value one of OR, PR, HR, beta), `exposure_scale` ("per 1 unit", "per SD", "Q4 vs Q1", or as written; if the estimate was fitted only in a restricted range or subset, say so, e.g. "per 1 unit, among TyG >= 9.05"), `estimate` (number), `ci_low` (number), `ci_high` (number), `p` (number or string as written), `source_table` (the table holding the estimate; its loc must be that table).
- `covariates_footnote` (field, value = list of covariates exactly as named in the footnote of the table that holds the main estimate).
- `covariates_methods` (field, value = list of covariates as named in the Methods text for that same model).
- `covariates_final` (plain list, no quote): the footnote list if a footnote list exists, otherwise the Methods list.
- `covariate_sources_differ` (boolean): true if the footnote list and the Methods list are not the same set of covariates (ignore only wording and ordering, not a missing or extra variable).
- `other_models` (list of objects `{"model_label","covariates","estimate","ci_low","ci_high","loc"}`): every other model in the same table for the same exposure and outcome, with its own covariate list as far as the paper states it.
- `inflection_point` (field: the threshold or inflection value if the paper reports one, else null).
- `ambiguities` (list of `{"point","reading_1","reading_2"}`).
- `notes` (free text, observations only).

## 2. Rule R1: which estimate is "the main estimate"

Do these steps in order.

1. Read the ABSTRACT. List every adjusted association estimate between the index and a kidney outcome that the abstract reports with a number.
2. Choose the outcome by this priority: CKD (including DKD or CKD in diabetes) first; if none, albuminuria; if none, low eGFR. Within the same outcome type, take the one stated first in the abstract.
3. Within that outcome, take the MOST FULLY ADJUSTED estimate shown in the abstract (the abstract often names the final model, or gives only one).
4. Prefer the continuous form of the index (per unit, per SD, log-transformed continuous) over quartiles, tertiles or categories. If the abstract gives only a categorical estimate for that outcome, use it; if it gives both, use the continuous one.
5. Find that same estimate in a results table. Take `model_label`, the estimate and CI, and the covariate footnote from THAT table. If the same number appears in several tables, use the one for the full sample, not a subgroup, sensitivity or stratified table.
6. If the abstract has no numeric estimate, or the abstract is the only text you have, apply the same logic to the first results table with the full-sample, most adjusted model, record what you did in `notes`, and add an ambiguity entry.
7. Put the outcome type chosen in `outcome_main` with a quote from the abstract or table row.

Do not choose an estimate because it is significant, large or interesting. Choose by this rule only.

## 3. Covariates: footnote versus Methods

- `covariates_footnote` comes from the footnote or legend of the table that holds the main estimate (words like "adjusted for", "Model 3 adjusted for age, sex, ..."). Copy the variable names as written. If models are nested ("Model 2 = Model 1 + ..."), expand the nesting into the full list for the chosen model in `value`, quote the footnote sentence that defines the chosen model, and say in `notes` that it was expanded.
- `covariates_methods` comes from the Methods (statistical analysis) description of that same model.
- Do not copy one into the other. Read both separately. If the paper has only one of them, leave the other as null and set `covariate_sources_differ` to false.
- If a covariate appears only as a group name ("socioeconomic status", "lifestyle factors", "comorbidities"), record it as written; do not expand it by guessing.
- A variable the sample was restricted on is not a covariate.
- Do not add the exposure, its components, or the outcome to the list unless the paper lists them as adjusted covariates.

## 4. Ambiguity readings

When the text leaves a specification point unclear in a way that would change what a reader would run (for example: which weight, how an outcome is defined when two definitions appear, whether "CKD" means eGFR below 60 only or also albuminuria, which eGFR equation, whether a covariate is categorical or continuous, whether a subpopulation is restricted), record it in `ambiguities`.

- Each entry: `point` (short name), `reading_1`, `reading_2`.
- `reading_1` is the MOST LITERAL reading of the words in the paper, even if you suspect the authors meant something else. `reading_2` is the next most plausible alternative.
- At most 2 readings per point, at most 3 points in total. Keep only the points that would most change an estimate. If a point has only one sensible reading, do not list it.
- Do not use ambiguities for typos or for things the paper simply does not report; use null fields for those.

## 5. Special cases

- Full text only partly available (abstract only, or tables missing): extract what can be found; set the missing fields to null; say so in `notes`. Never fill a table value from memory.
- Table values garbled by text extraction (columns run together): if you can still read the estimate and CI with certainty, quote the row exactly as it appears in the text (even if garbled); if not, set the estimate fields to null and say so in `notes`.
- The effect is reported as beta or on another scale: record effect_type as written and the scale as written.
- Several indices in one paper: the index follows from rule R1, not from the order of the objective sentence. The main index is the one whose estimate is picked by R1 (the first kidney-outcome estimate with a number in the abstract, in the order the abstract states them). Only when the abstract has no numeric estimate, take the index named first in the abstract objective (or title). Record the other candidate in `ambiguities`.
- Several NHANES cycle sets: list all cycles used in the main analysis.
- Multiple sources (for example NHANES plus a hospital cohort): extract the NHANES part only; mention this in `notes`.
- CI written as "(1.10-1.50)" or "1.10, 1.50": read the two bounds carefully; a hyphen is a separator, not a minus sign.

## 6. Before you answer

Check, silently:
1. Every `quote` can be found in the text exactly as you wrote it.
2. `estimate` lies between `ci_low` and `ci_high`, and all three belong to the same row of the same table as the quote.
3. The covariate footnote and the main estimate come from the same table (same table number in `loc`).
4. `outcome_main` follows rule R1.
5. The JSON is valid.

Then output the JSON object.

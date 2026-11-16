# I139 extraction output schema (shared by pass A and pass B)

One JSON object per paper, written to `<pass_dir>/<rec_id>.json`. Every key field is an object
`{"value": ..., "quote": "<verbatim text copied from the full text>", "loc": "<e.g. Methods 2.3 para 2 | Table 2 | Table 2 footnote | Abstract>"}`.
Use `"value": null, "quote": "", "loc": ""` when the paper does not report it. Never paraphrase inside `quote`; copy characters exactly (the checker matches after normalising whitespace and dashes only).

Top-level keys:
- `rec_id`, `doi` (plain strings)
- `index_name` (field), `index_formula` (field, value as written), `index_components` (field, value = list of component names, e.g. ["TG","FPG"])
- `cycles` (field, value = list like ["2011-2012","2013-2014"])
- `population` (field: age limits and main inclusion/exclusion as written)
- `n_analytic` (field, integer)
- `outcomes_in_abstract` (field, value = list of outcomes named in the abstract)
- `outcome_main` (field: the outcome of the main estimate, chosen by rule R1 below)
- `outcome_definition` (field), `egfr_equation` (field), `uacr_threshold` (field)
- `weight_variable` (field: name or description of the survey weight), `cycle_pooling` (field: how weights were combined across cycles)
- `main_estimate` (object): `model_label` (field, e.g. "Model 3"), `effect_type` (field: OR|PR|HR|beta), `exposure_scale` (field: "per 1 unit" | "per SD" | "Q4 vs Q1" | other as written), `estimate` (field, number), `ci_low` (field, number), `ci_high` (field, number), `p` (field, number or string as written), `source_table` (field: table holding the estimate; loc must be that table)
- `covariates_footnote` (field, value = list of covariates as written in the footnote of the table holding the main estimate)
- `covariates_methods` (field, value = list as written in Methods for that model)
- `covariates_final` (value = list; = footnote list if present, else methods list; plain list, no quote)
- `covariate_sources_differ` (boolean)
- `other_models` (list of {"model_label","covariates","estimate","ci_low","ci_high","loc"} for every other model in the same table)
- `inflection_point` (field: value, null if none)
- `ambiguities` (list of {"point","reading_1","reading_2"}; reading_1 = most literal = "first combination"; at most 2 readings per point, at most 3 points)
- `notes` (free text, observations only)

Rule R1 (prereg M1): the main estimate is the most fully adjusted estimate reported in the ABSTRACT, preferring the continuous index. If the abstract has several outcomes: CKD first, then albuminuria, then low eGFR; within the same type, the first one stated. Then locate that same estimate in the results table and take model label and covariates from that table.

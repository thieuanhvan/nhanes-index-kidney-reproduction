# Blind three-field audit: brief (written 15/11/2026 before any audit run)

You check data extracted from published papers against the papers' full text. You are an auditor; you do not know
why the data were extracted and you must not try to find out.

Working folder: <workdir>/idea-hunt-metaresearch-20261114

Read ONLY:
- `extraction/blind_audit/packet_blind_audit_i139.json` (the items assigned to you)
- the full text file named in each item (`fulltext_txt/<rec_id>.txt`; if a formula or table seems missing,
  `fulltext/<rec_id>.xml` may be read for that paper)

Do NOT read any other file in this folder or elsewhere (no extraction tables, labels, notes, verdicts, scripts,
results, prereg, scratchpad). Do not run models. Do not use the web.

Definition of the paper's "main estimate": the most fully adjusted estimate stated in the ABSTRACT, preferring the
continuous index over categories. If the abstract has several outcomes: CKD first, then albuminuria, then low eGFR;
within the same type, the first one stated. The value is then found in the body (usually a regression table).

For each item check three fields:
1. `model`: is this the label of the model (as named in the paper, e.g. "Model 3") that produced the main estimate?
2. `covariates`: is this the covariate list of that model as stated in the paper (table footnote or Methods)?
   Judge on content, not wording: synonyms and abbreviations are fine; a missing or extra variable is an error.
   Ignore the exposure itself.
3. `estimate`: are the OR (or effect estimate) and its 95% CI lower and upper bounds those of the main estimate?

Some fields may contain errors; most do not. Judge each field independently from the paper.

Output: write a TSV at the path given in your task with columns
`item, rec_id, field, verdict, correct_value, evidence_line`
- one row per item per field (3 rows per item)
- verdict = `correct` / `wrong` / `cannot_tell`
- correct_value: when wrong, what the paper says (for covariates: list the missing and extra variables as
  `missing: ...; extra: ...`); otherwise empty
- evidence_line: line number(s) in the full text file you relied on

Then report only: number of items done, counts of each verdict per field. Do not give opinions on the papers.

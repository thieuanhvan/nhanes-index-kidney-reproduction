# Reproducing published NHANES index and kidney disease studies: code and derived data

This repository holds the code, protocol, directed acyclic graphs (DAGs), sampling frame, derived per-paper results and run logs for the study "Reproducing published NHANES index and kidney disease studies from raw data: a pre-registered meta-research study of a stratified sample of 47 papers" (Van Thieu, University of Information Technology, Ho Chi Minh City, Vietnam). The protocol is registered on the Open Science Framework (https://osf.io/egjy3). Each folder keeps the name and layout it had during the analysis, so the relative paths inside the scripts work without changes.

## What is and is not included

Included: all analysis code, the registered protocol (`prereg_i139_v7_en.md`, byte-identical to the registered file; its hash is in the `.sha256` file), the DAGs, the registered sampling frame, the locked specification table, the per-paper results of every analysis stage, the prompts and briefs given to language-model agents, the reconciliation and amendment rules, and the run logs.

Not included:

- full texts, abstracts and verbatim sentences of the audited papers (copyright). Fields that held verbatim text are replaced by `[excerpt withheld]`, and table columns holding titles, quotes or evidence sentences are dropped;
- the raw extraction passes and the blind-audit packets, which consist largely of verbatim text;
- NHANES files, which are public at the CDC; `tools/download_nhanes_core.py` and `tools/nhanes_extra_i139.py` fetch them.

The audited papers are identified only by `rec_id`, a short hash of the record key. DOIs are withheld in every file except the registered frame (`frame/frame_i139.tsv`). The frame, the seed in the protocol and the sampling rule make the sample reproducible, so the papers can be identified by anyone who repeats the sampling; the paper does not name them.

Because DOIs and excerpts were removed, files whose hashes were recorded during the study (`extraction/spec_table_i139.sha256`, `extraction/labels_i139.sha256`, `extraction/m2_m3_rules_i139.sha256` and the audit locks) do not match the copies here; those hashes refer to the unredacted files kept by the author. The protocol, the frame and the DAG files are unredacted and match their hashes.

## Folders

| Folder | Contents |
|---|---|
| `tools/` | All scripts; `tools/m1_specs/<rec_id>.py` holds the per-paper reconstruction (first combination of readings of the locked specification) |
| `dag/` | DAG code (dagitty) and figures for the five DAG families, versions A (minimal) and B (extended), and the manifest |
| `frame/` | Registered sampling frame and its hash |
| `screening_full/` | Screening prompts, decisions of the two language-model passes (quote column dropped), disagreements and their reconciliation, final decisions, sampling queue |
| `extraction/` | Extraction prompts and briefs, reconciliation rules, the locked specification table, labels, results of M1 (reproduction), M2 (survey design), M3 (confounder specification) and M4 (inflection points), sensitivity analyses, overlap of participants, agreement statistics, audit results, logs |
| `tests/` | Unit tests for the index reconstruction |

## Run order

1. Frame and screening: `build_frame_i139.py`, `fetch_screen_i139.py`, `build_packet_i139.py`; two screening passes with `screening_full/prompt_screen_A.md` and `prompt_screen_B.md`; `screen_agreement.py`.
2. Full texts and extraction: `fetch_fulltext_i139.py`, `to_text_i139.py`; two extraction passes with `extraction/prompt_extract_A.md` and `prompt_extract_B.md`; `check_extraction_i139.py`, `compare_extraction_i139.py`, `kappa_ac1_i139.py`; reconciliation under `extraction/reconcile_brief_main.md`.
3. Specification and labels: `spec_table_i139.py`, `concepts_i139.py`, `labels_i139.py`.
4. NHANES: `download_nhanes_core.py`, `nhanes_extra_i139.py`, `nhanes_core.py`, `nhanes_indices.py`, `diet_indices_i139.py`, `recon_i139.py`.
5. Analyses: `m1_run_i139.py` (with `m1_lib_i139.py` and `m1_specs/`), `m2_run_i139.py`, `m3_run_i139.py`, `m3_summary_i139.py`, `m3_simpc_i139.py`, `m3_desc_i139.py`, `m3_outsens_i139.py`, `m3_twoorder_i139.py`, `m4_run_i139.py`, `jaccard_i139.py`.
6. Verification: `inject_errors_i139.py`, `audit_openai_i139.py`, `kd3_variance_check.py`.
7. Analyses added after the main results were seen (labelled post hoc in the paper): `m1_anycomb_i139.py`, `m1_fail_causes_i139.py`, `sens_sample_i139.py`, `r3_m3m4_i139.py`, `r3_weighted_primary_i139.py`; competing-work scan `kd5_presubmission_scan_i139.py`.

Python 3 with numpy, pandas, scipy and statsmodels. Steps 1 and 2 call language-model agents and web services, so they cannot be repeated exactly; their outputs from this study are in `screening_full/` and `extraction/`. Steps 3 to 7 are deterministic given the specification table and the NHANES files.

`tools/build_human_audit_i139*.py` built the form for a human re-extraction of 32 papers that was planned but not done (deviation 2 in Supplementary Table S1 of the paper).

## Licence

Code: MIT License. Derived data, DAGs and documents: CC BY 4.0. See `LICENSE`.

# CUP survival surveillance analysis code

Statistical code for registry-defined glandular cancer of unknown primary (CUP)
and distant-stage known-primary lung, pancreatic, colorectal and stomach cancers.
The analysis examines first-to-last-era absolute survival-gap changes, early
mortality, conditional survival and sensitivity analyses using SEER 17 data.

**This repository contains code only, plus configuration and operating instructions.**
It contains no manuscript, patient-level records, aggregate results, figures,
model objects, bootstrap draws, SEER export files or personal contact details.
Obtain your own authorized data access; this repository does not provide access
to SEER data. See <https://seer.cancer.gov/data/access.html>.

## Software and design

- Data extraction: SEER*Stat 9.0.43.0.
- Statistical runtime used: R 4.6.1; data.table 1.18.4 and jsonlite 2.0.0.
- R base/recommended packages: stats, splines, parallel, survival, grid.
- Optional figure modules additionally use ggplot2, patchwork, ragg and svglite.
- The incidence verifier and repository audit require only the Python standard library.
- Scripts do not install dependencies, download data, or upload generated files.

Membership odds standardize each source group to the measured CUP reference
distribution. Weighted Kaplan–Meier curves provide 90-day and 1-year survival,
conditional survival and restricted mean survival time. Shared patient bootstrap
multiplicities preserve overlapping source/reference roles; seed for draw `b`
is `20261025 + b`. The analysis uses 2,000 draws, fixed initial spline bases,
nominal secondary intervals and four-comparison Bonferroni primary intervals.
Morphology definitions and primary contrasts were designated after preliminary
analyses; additional sensitivity analyses are exploratory. Code publication is
not a claim of prospective registration, causal identification or independent validation.

## Folder map

- `R/method_reproducibility/`: cohort preparation, fitting kernel, reference analyses.
- `R/validation_stage3/`: broad/strict-definition bootstrap analyses.
- `R/final_audit/`: four-era and annual broad/strict analyses and year mixtures.
- `R/histology_revision/`: glandular definitions, primary models, paired contrasts,
  year mixtures, targeted sensitivity analyses and interval checks.
- `R/pandemic_sensitivity/`: contrasts excluding pandemic-era diagnoses.
- `R/clinical_readability/`: absolute outcomes and cohort composition.
- `R/weight_capping/`: source-odds winsorization with balance reassessment.
- `R/early_landmark/`: early outcomes, survivor case mix and landmark standardization.
- `R/tissue_inference/`: histologic-confirmation restriction and paired inference.
- `R/precision_polish/`: early-outcome figure-generation code.
- `config/histology_definitions.json`: code lists and statistical parameters, not data.
- `python/verify_incidence.py`: observed age-adjusted Rate export checks.
- `tests/smoke.R`: in-memory synthetic computational checks.

Names of the analysis folders and some output basenames are retained for
traceability to the validated implementation. They do not contain archived
analysis outputs. Code was made portable by replacing machine paths, relocating
script/configuration references, and making worker count configurable. Statistical
formulae, morphology code lists, seed rules and inferential thresholds were retained.

## Input preparation

See [INPUT_SCHEMA.md](INPUT_SCHEMA.md). Keep your authorized exports **outside**
this repository. Set `CUP_EXPORT_DIR` to the folder containing the five tab-delimited
Case Listing exports. `CUP_PROJECT_ROOT` defaults to the repository root.

Example in an R session, replacing the placeholder with your private input folder:

```r
Sys.setenv(CUP_EXPORT_DIR = "<PRIVATE_EXPORT_DIRECTORY>", CUP_WORKERS = "2")
system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "run_analysis.R", "prepare"))
```

Do not put credentials in a script, commit an environment file, or force-add data.
The allowlist `.gitignore` excludes generated directories and all data formats by default.

## Run order

From the repository root:

```sh
Rscript --vanilla tests/smoke.R
Rscript --vanilla run_analysis.R list
Rscript --vanilla run_analysis.R all
```

`all` is computationally intensive: several analyses refit weights in each of
2,000 draws. `CUP_WORKERS` controls PSOCK workers (default two). Each stage runs
in a fresh R process. Intermediate checkpoints and results are written locally
to ignored `work/` and `outputs/` directories. Patient-level intermediate objects
must remain private. Do not change data, code, seed or design while resuming a checkpoint.

For staged execution, choose names printed by `list`. Dependencies follow the
listed order. `full_record_check`, `broad_fit`, `broad_summary`, `glandular_fit`
and `glandular_summary` supply the earlier-definition fits needed by paired
morphology comparisons. Later reference, tissue and early-outcome analyses use
the generated common bootstrap indexing. `reference_reproduction` is required
before `early_outcomes`. Figure modules are optional; no precomputed figure is supplied.

## Incidence export

In SEER*Stat, select age 20+, diagnoses 2004–2022, malignant C80.9, sequence
00/01, microscopic confirmation 1/2/4 and the 53 retained morphology codes.
Do not impose survival eligibility on the incidence cohort. Match the population
and standard-population age/year filters. Use direct age adjustment to the
2000 U.S. standard population, adult strata renormalized to sum to one, rates
per 100,000, and Tiwari-modified 95% intervals. For a delay-factor database,
select **Delay site = All sites** and **Delay race = All Races** to prevent
duplicate representations. Export ordinary observed Rate/SE/CI/Count/Population
fields, rather than delay-adjusted estimates for this restricted case definition.

```sh
python python/verify_incidence.py <PRIVATE_AGGREGATE_RATE_EXPORT>
```

The verifier checks the exported estimates; it cannot reconstruct age-specific
rates from an aggregate table. It does not fit Joinpoint, APC or AAPC models.
Generated incidence files remain ignored. Official age-adjustment documentation:
<https://seer.cancer.gov/seerstat/tutorials/aarates/definition.html>.

## Verification limits and release hygiene

The public synthetic test does not reproduce study estimates without the licensed
inputs. Before publication, the fitting kernel was checked against the retained
implementation and cohort point estimates using local authorized inputs; those
inputs and verification outputs are not published. The full set of long-running
bootstrap analyses was not rerun solely to prepare this portable code release.

Before any later push, stage only intended files and run:

```sh
python tools/audit_code_only.py
git diff --cached --name-only
```

The scanner is an additional check, not a guarantee that arbitrary future files
are safe to disclose. Never bypass the allowlist with `git add --force` for data.

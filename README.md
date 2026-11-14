# CUP survival analysis

R code for the final analysis of registry-defined glandular cancer of unknown
primary (CUP) and distant-stage lung, pancreatic, colorectal and stomach cancers
in SEER Research Data, 17 registries, November 2025 submission. The study covers
diagnoses in 2004–2022 and follow-up through December 31, 2023.

## Analyses

The pipeline constructs the eligible cohort, fits membership-odds weights,
checks balance and effective sample size, and estimates weighted Kaplan–Meier
survival, restricted mean survival time and conditional survival. Each cancer
and diagnosis era is standardized to a common pooled CUP reference.
The primary estimates are the four changes in CUP-minus-benchmark 1-year
survival gaps between 2004–2007 and 2020–2022.

The analysis includes morphology and histologic-confirmation restrictions,
annual and pre-pandemic comparisons, an early-era CUP reference, weight
winsorization, survivor case mix and landmark standardization. Clinical analyses
provide a symmetric decomposition of 1-year survival changes, age strata
(<65, 65–74, ≥75), paired between-stratum contrasts, and a 60-day survival
boundary. CSS uses the SEER cause-specific death classification under two rules
for unknown cause of death.

Inferential modules use 2,000 shared patient bootstrap draws, stratified
by cancer and era, with seed 20261025 + b for draw b. Membership models
are refitted in each draw with fixed initial spline bases. Primary intervals
use four-comparison Bonferroni percentiles. Exploratory clinical and CSS
intervals are nominal 95% percentile intervals. Primary estimands, definitions
and some analytic choices were designated after initial data inspection.
The clinical, age, boundary and CSS analyses are exploratory.

## Software

The analysis used R 4.6.1, data.table 1.18.4 and jsonlite 2.0.0, together with
base/recommended packages including splines, parallel and survival.
Figure scripts use ggplot2, patchwork, ragg and svglite.
Python scripts use the standard library. SEER*Stat 9.0.43.0 was used for extraction.
Use existing installations of these dependencies.

## Inputs and execution

Obtain authorized SEER data and prepare exports described in
[INPUT_SCHEMA.md](INPUT_SCHEMA.md). Store them in private directories.
Set CUP_EXPORT_DIR to the private cohort export folder, CUP_CSS_EXPORT_DIR
to the private death-classification export folder, and CUP_WORKERS to the
desired PSOCK worker count (default 2).

Run the following commands from the repository root, in an environment
with those variables set.

    Rscript --vanilla tests/smoke.R
    Rscript --vanilla run_analysis.R list
    Rscript --vanilla run_analysis.R all
    python tools/audit_code_only.py

CUP_PROJECT_ROOT defaults to the repository root. Stages run in separate
R processes in dependency order. For staged execution, use names printed
by list and complete preceding stages first. The full run is computationally
intensive. Checkpoints and generated outputs are local, ignored artifacts.
Resume a checkpoint with the same cohort, code, configuration and seed.

## Code organization

| Module | Purpose |
|---|---|
| R/standardization | Cohort preparation, fitting kernel and reference standardization |
| R/cohort_checks | Independent cohort and fixed-target point checks |
| R/definition_models | Broad and strict morphology models and paired draws |
| R/era_models | Four-era and annual models and year mixtures |
| R/glandular_models | Final glandular definition, primary contrasts and definition sensitivities |
| R/pandemic_sensitivity | Pre-pandemic diagnosis contrasts |
| R/absolute_outcomes | Absolute probabilities and cohort composition |
| R/weight_capping | Winsorized weights with balance reassessment |
| R/early_landmark | Early outcomes, survivor case mix and landmark standardization |
| R/tissue_inference | Histologic-confirmation restriction and paired contrasts |
| R/clinical_analyses | Survival decomposition, age strata and the 60-day boundary |
| R/cause_specific | Death-classification coverage and CSS sensitivity |
| R/figures | Main survival-gap and early-outcome figures |
| config/histology_definitions.json | Morphology codes and statistical parameters |
| python/verify_incidence.py | Checks of an authorized aggregate Rate export |

## Interpretation

The decomposition expresses the identity S(365.25) = S(90) ×
S(365.25)/S(90), using unrounded estimates and paired draws. Its components
describe the mathematical contribution of each factor. Age-stratum comparisons
use separate pooled CUP reference populations. Support requires source n ≥100,
ESS ≥200 and maximum absolute model-matrix SMD <0.10. The ≥75 stomach
model for 2020–2022 fails support and its estimates require caution.

CSS treats attributed cancer deaths as events and other-cause deaths as censored.
Unknown-cause deaths are either censored or counted as cancer deaths.
Coverage uses all-follow-up deaths as the denominator and also reports first-year
counts. Death attribution and censoring assumptions remain limitations.
The CSS module independently checks point estimates against survival::survfit
and matches its OS bootstrap estimates to the primary draws.

For incidence, use a separately extracted age-adjusted SEER Rate session
with matched population denominators. The restricted definition uses age ≥20,
C80.9, malignant behavior, sequence 00/01, microscopic confirmation 1/2/4
and the 53 glandular codes. See the input schema and incidence verifier.

## Data access

This repository publishes analysis code, configuration and operating instructions.
SEER access is governed by its data-use agreement.
See <https://seer.cancer.gov/data/access.html>.
Individual records, exports, model objects, bootstrap draws, aggregate results,
figures and manuscript files remain in private local directories.

# Private inputs and cohort construction

No data are included. Users must obtain the SEER Research Data under the
applicable data-use agreement and create their own SEER*Stat exports.

## Five Case Listing inputs

The preparation script expects tab-delimited exports with a header, in this order:

1. `S2_CUP_Adeno_2004_2022_labels.txt`
2. `S2_Lung_Adeno_Distant_2004_2022_labels.txt`
3. `S2_Pancreas_Adeno_Distant_2004_2022_labels.txt`
4. `S2_Colorectal_Adeno_Distant_2004_2022_labels.txt`
5. `S2_Stomach_Adeno_Distant_2004_2022_labels.txt`

Use the same 26 ordered columns in all five files. The script imports character
labels and renames them to the internal names below:

| Position | SEER field | Internal name |
|---|---|---|
| 1 | Patient ID | id |
| 2 | Record number recode | record |
| 3 | Sequence number | sequence |
| 4 | Year of diagnosis | year |
| 5 | Age recode with single ages and 85+ | age_label |
| 6 | Age recode with <1 year olds and 90+ | age_band |
| 7 | Sex | sex |
| 8 | Race and origin recode | race |
| 9 | Marital status at diagnosis | marital |
| 10 | Primary Site - labeled | site |
| 11 | Histologic Type ICD-O-3 | hist |
| 12 | Behavior code ICD-O-3 | behavior |
| 13 | Diagnostic Confirmation | dc |
| 14 | Type of Reporting Source | source |
| 15 | Survival months | months |
| 16 | Survival months flag | flag |
| 17 | Survival Days | days |
| 18 | Vital status recode (study cutoff used) | vital |
| 19 | Year of follow-up recode | followyear |
| 20 | Median household income inflation-adjusted to 2024 | income |
| 21 | Rural-Urban Continuum Code | rural |
| 22 | Combined Summary Stage with Expanded Regional Codes (2004+) | stage |
| 23 | SEER Combined Mets at DX-bone (2010+) | bone |
| 24 | SEER Combined Mets at DX-brain (2010+) | brain |
| 25 | SEER Combined Mets at DX-liver (2010+) | liver |
| 26 | SEER Combined Mets at DX-lung (2010+) | lung |

Export labels, not numeric recodes, except the numeric morphology field. Exact
labels can vary between SEER releases; inspect them locally before using a
different release. The analysis used the SEER 17 November 2025 submission.

## Initial export filters and analytic eligibility

Initial Case Listing exports: diagnoses 2004–2022, age at least 20, malignant
behavior, and morphology ranges 8140–8239, 8250–8384, 8440–8490, 8550–8551,
8570–8574 and 8576. The broad extraction is needed to reproduce paired
definition sensitivities; do not apply only the final whitelist at extraction.

CUP uses C80.9. Known-primary controls use C34.*, C25.*, C18.* with C19.9/C20.9,
and C16.*, respectively, and require distant disease in the combined summary-stage
variable. Appendiceal primaries are retained in the colorectal extraction and
excluded in a sensitivity analysis. CUP is not filtered by distant summary stage.

Preparation keeps sequence 00/01, excludes autopsy-only/death-certificate-only
records, requires microscopic confirmation 1/2/4, and valid survival/vital status.
The final glandular and sensitivity definitions are applied by the analysis
modules using `config/histology_definitions.json`.

Age is top-coded at 85 in the single-age field; a separate indicator identifies
90+ from the grouped field. Death is the event, zero-day survival is set to
0.5 day, and the 1-year horizon is 365.25 days. Unknown categories are retained.
The initial broad eligible file is the shared bootstrap sampling frame.

## Aggregate Rate input

The optional incidence verifier expects tab-separated quoted headers:
`Year of diagnosis`, `Age-Adjusted Rate`, `Standard Error`,
`Lower Confidence Interval`, `Upper Confidence Interval`, `Count`, `Population`,
`Std Pop`. Export a total row and annual rows. Only 2004–2022 enter the annual
series; out-of-range zero-population rows are omitted. This is a separate
population-denominator analysis, not a rate calculated from survival-eligible cases.

## Cause-specific death classification inputs

Set CUP_CSS_EXPORT_DIR to five tab-delimited files with headers exported from
the same SEER release and case definitions as the case listing inputs.

| Group | Filename |
|---|---|
| CUP | S2_CUP_Adeno_2004_2022_labels.txt |
| Lung | S2_Lung_Adeno_Distant_2004_2022_CSS.txt |
| Pancreas | S2_Pancreas_Adeno_Distant_2004_2022_CSS.txt |
| Colorectal | S2_Colorectal_Adeno_Distant_2004_2022_CSS.txt |
| Stomach | S2_Stomach_Adeno_Distant_2004_2022_CSS.txt |

Required headers are Patient ID, Record number recode, and
SEER cause-specific death classification. Other exported columns are permitted.
The CUP file can share its basename with the cohort export and reside
in the CSS input directory.

Patient ID and record number are imported as character values and uniquely
identify the matching record. The module requires complete classification
coverage of the final glandular cohort, with no duplicate input keys.
Alive patients must have classification 0.

Accepted classification labels and numeric values are 0 (Alive or dead of
other cause), 1 (Dead attributable to this cancer dx), 8 (Dead with
missing/unknown COD), and 9 (N/A not seq 0-59). The script uses the exact SEER
labels in its coding map. Code 9 and missing/unrecognized values fail the
analytic-cohort check. Use the published classification directly.
See https://seer.cancer.gov/causespecific/ for classification rules.

CSS censors other-cause deaths at their recorded death times.
For code 8, separate analyses censor at death or treat death as cancer death.
First-year estimates use 365.25 days and the same baseline weights and paired
patient draws as OS.

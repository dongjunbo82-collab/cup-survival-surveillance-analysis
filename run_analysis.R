# Run from the repository root. Data must be obtained separately under SEER terms.
options(stringsAsFactors=FALSE)
root <- normalizePath(Sys.getenv('CUP_PROJECT_ROOT', unset=getwd()), winslash='/')
setwd(root); Sys.setenv(CUP_PROJECT_ROOT=root)
steps <- list(
  prepare='R/standardization/prepare_eligible.R',
  full_record_check='R/cohort_checks/fixed_target_era_check.R',
  broad_fit='R/definition_models/fit_definition_models.R',
  broad_summary='R/definition_models/summarize_definition_models.R',
  four_era_fit='R/era_models/fit_era_models.R',
  four_era_summary='R/era_models/summarize_era_models.R',
  glandular_fit='R/glandular_models/fit_glandular_models.R',
  glandular_summary='R/glandular_models/summarize_glandular_models.R',
  definition_sensitivities='R/glandular_models/definition_sensitivities.R',
  interval_precision='R/glandular_models/check_interval_precision.R',
  pandemic='R/pandemic_sensitivity/pandemic_sensitivity.R',
  absolute_outcomes='R/absolute_outcomes/absolute_outcomes.R',
  early_reference='R/standardization/early_reference_sensitivity.R',
  reference_reproduction='R/standardization/reproduce_first_last.R',
  weights='R/weight_capping/weight_sensitivity.R',
  landmark='R/early_landmark/early_landmark_diagnostics.R',
  early_outcomes='R/early_landmark/early_era_summary.R',
  histologic_confirmation='R/tissue_inference/tissue_inference.R',
  clinical_analyses='R/clinical_analyses/clinical_analyses.R',
  cause_specific='R/cause_specific/cause_specific.R',
  survival_figures='R/glandular_models/make_figures.R',
  annual_composition='R/standardization/plot_annual_composition.R',
  early_outcomes_figure='R/figures/plot_early_outcomes.R',
  survival_gaps_figure='R/figures/plot_survival_gaps.R'
)
args <- commandArgs(TRUE)
if (!length(args) || identical(args,'list')) {
  cat(paste(names(steps),unlist(steps),sep=': ',collapse='\n'),'\n')
  quit(save='no')
}
requested <- if(identical(args,'all')) names(steps) else args
stopifnot(all(requested %in% names(steps)))
for (step in requested) {
  message('Running ',step)
  status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(steps[[step]])))
  if (status != 0L) stop('Step failed: ',step,'; inspect local outputs before continuing.')
}

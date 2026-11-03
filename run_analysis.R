# Run from the repository root. Data must be obtained separately under SEER terms.
options(stringsAsFactors=FALSE)
root <- normalizePath(Sys.getenv('CUP_PROJECT_ROOT', unset=getwd()), winslash='/')
setwd(root); Sys.setenv(CUP_PROJECT_ROOT=root)
steps <- list(
  prepare='R/method_reproducibility/prepare_eligible.R',
  full_record_check='R/stage2_preanalysis/fixed_target_era_check.R',
  broad_fit='R/validation_stage3/run_validation.R',
  broad_summary='R/validation_stage3/summarize_validation.R',
  four_era_fit='R/final_audit/run_final_audit.R',
  four_era_summary='R/final_audit/summarize_final_audit.R',
  glandular_fit='R/histology_revision/run_revision.R',
  glandular_summary='R/histology_revision/summarize_revision.R',
  definition_sensitivities='R/histology_revision/revision_sensitivities.R',
  interval_precision='R/histology_revision/check_interval_precision.R',
  pandemic='R/pandemic_sensitivity/pandemic_sensitivity.R',
  absolute_outcomes='R/clinical_readability/clinical_readability.R',
  early_reference='R/method_reproducibility/early_reference_sensitivity.R',
  reference_reproduction='R/method_reproducibility/reproduce_first_last.R',
  weights='R/weight_capping/weight_sensitivity.R',
  landmark='R/early_landmark/early_landmark_diagnostics.R',
  early_outcomes='R/early_landmark/early_era_summary.R',
  histologic_confirmation='R/tissue_inference/tissue_inference.R',
  survival_figures='R/histology_revision/make_figures.R',
  annual_composition='R/method_reproducibility/plot_annual_composition.R',
  early_outcomes_figure='R/precision_polish/plot_S4_early_outcomes.R'
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

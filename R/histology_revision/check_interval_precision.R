# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
library(data.table)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'))
out<-file.path(root,'outputs/histology_revision')
A<-readRDS(file.path(root,'work/histology_review/summary_arrays.rds'))$A
z<-rbindlist(lapply(c('Lung','Pancreas','Colorectal','Stomach'),function(g){
 ls<-paste('Era','Reviewed',c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__')
 v<-A[,ls[1],'OS365']-A[,ls[2],'OS365']-A[,ls[3],'OS365']+A[,ls[4],'OS365']
 a<-quantile(v[1:1000],c(.00625,.99375));b<-quantile(v[1001:2000],c(.00625,.99375))
 data.table(control=g,first_half_lo=a[1],first_half_hi=a[2],second_half_lo=b[1],second_half_hi=b[2],maximum_endpoint_difference_pp=100*max(abs(a-b)),both_halves_exclude_zero=all(c(a,b)<0))
}))
fwrite(z,file.path(out,'32_primary_interval_precision.csv'),bom=TRUE)
stopifnot(all(z$both_halves_exclude_zero));print(z)

# Generated artifacts are written to local output directories.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/glandular_models', 'work/landmark',
            'outputs/standardization', 'outputs/cohort_checks',
            'outputs/definition_models', 'outputs/era_models', 'outputs/glandular_models',
            'outputs/pandemic_sensitivity', 'outputs/absolute_outcomes',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
# Rebuild both reference analyses from the authorized eligible file; using shared patient bootstrap sampling.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(splines);library(parallel);library(jsonlite)})
setDTthreads(1);root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'))
out<-file.path(root,'outputs/standardization');e<-readRDS(file.path(root,'work/eligible_cohort.rds'))
source(file.path(root,'R/standardization/analysis_functions.R'),encoding='UTF-8');cfg<-fromJSON(file.path(root,'config/histology_definitions.json'))
keep<-e$hist%in%cfg$glandular_codes;groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-c('2004-2007','2020-2022')
models<-list();points<-list();diag<-list()
for(ref in c('Pooled','Early')){
 target<-which(keep&e$group=='CUP'&(ref=='Pooled'|as.character(e$era)=='2004-2007'))
 for(er in eras)for(g in groups){lab<-paste(ref,er,g,sep='__');m<-build_model(target,which(keep&e$group==g&as.character(e$era)==er),lab,ref,er,g);f<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(f$error));models[[lab]]<-m;points[[lab]]<-data.table(reference=ref,era=er,group=g,label=lab,measure=measures,estimate=f$values);diag[[lab]]<-cbind(data.table(reference=ref,era=er,group=g),as.data.table(as.list(f$diagnostic)))}
}
pt<-rbindlist(points);fwrite(rbindlist(diag),file.path(out,'58_standalone_reference_diagnostics.csv'))
repfit<-function(b){set.seed(20261025+b);counts<-integer(N);for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx));v<-matrix(NA_real_,length(models),6);for(i in seq_along(models)){z<-fit_model(models[[i]],counts);if(nzchar(z$error))stop(z$error);v[i,]<-z$values};v}
B<-2000L;cp<-file.path(root,'work/landmark/standalone_reproduction.rds');sigfile<-paste0(cp,'.signature.rds')
dir.create(dirname(cp),recursive=TRUE,showWarnings=FALSE)
sig<-list(input=unname(tools::md5sum(file.path(root,'work/eligible_cohort.rds'))),functions=unname(tools::md5sum(file.path(root,'R/standardization/analysis_functions.R'))),config=unname(tools::md5sum(file.path(root,'config/histology_definitions.json'))),labels=names(models),B=B,seed=20261025)
if(file.exists(sigfile))stopifnot(identical(sig,readRDS(sigfile)))else saveRDS(sig,sigfile)
boot<-if(file.exists(cp))readRDS(cp)else list();stopifnot(length(boot)<=B)
cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS','2')),outfile=file.path(root,'work/landmark/standalone_workers.log'));clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL});clusterExport(cl,c('N','strata','models','repfit','fit_model','km_cached'),envir=environment())
for(start in if(length(boot)<B)seq(length(boot)+1L,B,by=100)else integer()) {boot<-c(boot,parLapplyLB(cl,start:min(B,start+99L),repfit));saveRDS(boot,cp);cat('Completed',length(boot),'/',B,'\n');flush.console()};stopCluster(cl)
stopifnot(length(boot)==B,all(vapply(boot,function(v)all(is.finite(v)),logical(1))))
sign<-c(1,-1,-1,1);res<-list()
for(ref in c('Pooled','Early'))for(g in groups[-1]){
 labs<-paste(ref,c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__');ii<-match(labs,names(models));v<-vapply(boot,function(b)sum(b[ii,2]*sign),numeric(1));p<-sum(vapply(labs,function(l)pt[label==l&measure=='OS365',estimate],numeric(1))*sign)
 res[[length(res)+1]]<-data.table(reference=ref,control=g,estimate_pp=100*p,lo_family_pp=100*quantile(v,.00625),hi_family_pp=100*quantile(v,.99375))
}
res<-rbindlist(res);fwrite(res,file.path(out,'57_standalone_reproduced_changes.csv'))
# Compare independent reference estimates with the common model summaries when available.
tests<-list(data.table(check='2000 bootstrap samples; 20 finite models in each',error=0,passed=TRUE))
f<-file.path(out,'48_early_reference_primary_changes_CI.csv')
if(file.exists(f)){z<-merge(res[reference=='Early'],fread(f),by='control',suffixes=c('_comparison','_reference'));err<-max(abs(as.matrix(z[,.(estimate_pp_comparison,lo_family_pp_comparison,hi_family_pp_comparison)])-as.matrix(z[,.(estimate_pp_reference,lo_family_pp_reference,hi_family_pp_reference)])));stopifnot(err<1e-10);tests[[length(tests)+1]]<-data.table(check='Standalone early-reference results match shared-bootstrap estimates and CIs',error=err,passed=TRUE)}
f<-file.path(root,'outputs/glandular_models/17_gap_changes_CI.csv')
if(file.exists(f)){z<-fread(f)[definition=='Glandular'&era=='2020-2022'&measure=='OS365',.(control,estimate_pp_reference=estimate*100,lo_family_pp_reference=lo_family*100,hi_family_pp_reference=hi_family*100)];z<-merge(res[reference=='Pooled'],z,by='control');err<-max(abs(as.matrix(z[,.(estimate_pp,lo_family_pp,hi_family_pp)])-as.matrix(z[,.(estimate_pp_reference,lo_family_pp_reference,hi_family_pp_reference)])));stopifnot(err<1e-10);tests[[length(tests)+1]]<-data.table(check='Standalone pooled-reference results reproduce primary estimates and CIs',error=err,passed=TRUE)}
fwrite(rbindlist(tests),file.path(out,'59_standalone_reproduction_QA.csv'));print(res)

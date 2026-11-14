# Generated artifacts are written to local output directories.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/glandular_models', 'work/landmark',
            'outputs/standardization', 'outputs/cohort_checks',
            'outputs/definition_models', 'outputs/era_models', 'outputs/glandular_models',
            'outputs/pandemic_sensitivity', 'outputs/absolute_outcomes',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(splines);library(parallel);library(jsonlite)})
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'))
e<-readRDS(file.path(root,'work/eligible_cohort.rds'))
source(file.path(root,'R/standardization/analysis_functions.R'),encoding='UTF-8');setDTthreads(1)
dest<-file.path(root,'outputs/standardization');dir.create(dest,showWarnings=FALSE)
cfg<-fromJSON(file.path(root,'config/histology_definitions.json'));keep<-e$hist%in%cfg$glandular_codes
target<-which(keep&e$group=='CUP'&as.character(e$era)=='2004-2007')
groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-c('2004-2007','2020-2022')
models<-list();points<-list();diags<-list();balances<-list();knots<-list()
for(er in eras)for(g in groups){
 source<-which(keep&e$group==g&as.character(e$era)==er);lab<-paste(er,g,sep='__')
 m<-build_model(target,source,lab,'EarlyReference',er,g);z<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(z$error))
 models[[lab]]<-m;points[[lab]]<-data.table(label=lab,era=er,group=g,measure=measures,estimate=z$values)
 w<-z$odds;diags[[lab]]<-cbind(data.table(label=lab,era=er,group=g),as.data.table(as.list(c(z$diagnostic,weight_min=min(w),weight_p99=unname(quantile(w,.99)),weight_mean=mean(w),top1pct_share=sum(sort(w,decreasing=TRUE)[seq_len(ceiling(length(w)*.01))])/sum(w)))))
 balances[[lab]]<-cbind(data.table(label=lab),z$smd)
 a<-ns(c(e$age[target],e$age[source]),df=4)
 knots[[lab]]<-data.table(label=lab,internal_knots=paste(attr(a,'knots'),collapse=';'),boundary_knots=paste(attr(a,'Boundary.knots'),collapse=';'),reference_source_overlap=length(intersect(target,source)))
}
pt<-rbindlist(points);dg<-rbindlist(diags);dg[,support_pass:=n_source>=100&ESS>=200&max_SMD<.1]
fwrite(pt,file.path(dest,'43_early_reference_point_estimates.csv'));fwrite(dg,file.path(dest,'44_early_reference_diagnostics.csv'));fwrite(rbindlist(balances),file.path(dest,'45_early_reference_balance.csv'));fwrite(rbindlist(knots),file.path(dest,'46_early_reference_spline_settings.csv'))
saveRDS(models,file.path(root,'work/landmark/early_reference_models.rds'))
# Verify reused fitting algorithm against the glandular primary model.
oldm<-readRDS(file.path(root,'work/glandular_models/models.rds'))[['Era__Glandular__2020-2022__Lung']]
oldp<-fread(file.path(root,'outputs/glandular_models/04_point_estimates.csv'))[label=='Era__Glandular__2020-2022__Lung',estimate]
stopifnot(max(abs(fit_model(oldm)$values-oldp))<1e-12)
repfit<-function(b){set.seed(20261025+b);counts<-integer(N);for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx));v<-matrix(NA_real_,length(models),6);err<-character(length(models));for(i in seq_along(models)){z<-tryCatch(fit_model(models[[i]],counts),error=function(ex)list(values=rep(NA_real_,6),error=conditionMessage(ex)));v[i,]<-z$values;err[i]<-z$error};list(b=b,v=v,error=err)}
B<-2000L;cp<-file.path(root,'work/landmark/early_reference_bootstrap.rds');boot<-if(file.exists(cp))readRDS(cp)else list()
cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS','2')),outfile=file.path(root,'work/landmark/early_reference_workers.log'));clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL});clusterExport(cl,c('N','strata','models','repfit','fit_model','km_cached'),envir=environment())
for(start in if(length(boot)<B)seq(length(boot)+1L,B,by=100)else integer()){
 boot<-c(boot,parLapplyLB(cl,start:min(B,start+99L),repfit));saveRDS(boot,cp);cat('Completed',length(boot),'/',B,format(Sys.time()),'\n');flush.console()
}
stopCluster(cl)
failure<-rbindlist(lapply(boot,function(z)data.table(replicate=z$b,label=names(models),error=z$error)));fwrite(failure[nzchar(error)],file.path(dest,'47_bootstrap_failures.csv'));stopifnot(!any(nzchar(failure$error)))
AA<-array(NA_real_,c(B,length(models),6),dimnames=list(NULL,names(models),measures));for(b in 1:B)AA[b,,]<-boot[[b]]$v
oldA<-readRDS(file.path(root,'work/glandular_models/summary_arrays.rds'))$A
oldch<-fread(file.path(root,'outputs/glandular_models/17_gap_changes_CI.csv'))[definition=='Glandular'&era=='2020-2022'&measure=='OS365']
res<-list();abs<-list()
for(l in names(models))for(k in seq_along(measures)){vv<-AA[,l,k];abs[[length(abs)+1]]<-data.table(label=l,measure=measures[k],estimate=pt[label==l&measure==measures[k],estimate],lo=quantile(vv,.025),hi=quantile(vv,.975))}
for(g in groups[-1]){
 labs<-paste(c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__');sgn<-c(1,-1,-1,1)
 p<-sum(vapply(labs,function(l)pt[label==l&measure=='OS365',estimate],numeric(1))*sgn)
 v<-as.numeric(AA[,labs,'OS365']%*%sgn);ol<-paste('Era','Glandular',c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__');ov<-as.numeric(oldA[,ol,'OS365']%*%sgn)
 res[[g]]<-data.table(control=g,estimate_pp=100*p,lo_pp=100*quantile(v,.025),hi_pp=100*quantile(v,.975),lo_family_pp=100*quantile(v,.00625),hi_family_pp=100*quantile(v,.99375),paired_shift_pp=100*(p-oldch[control==g,estimate]),paired_shift_lo_pp=100*quantile(v-ov,.025),paired_shift_hi_pp=100*quantile(v-ov,.975),support_pass=all(dg[label%in%labs,support_pass]))
}
res<-rbindlist(res);fwrite(res,file.path(dest,'48_early_reference_primary_changes_CI.csv'));fwrite(rbindlist(abs),file.path(dest,'49_early_reference_absolute_CI.csv'))
cup<-e[keep&group=='CUP']
annual<-cup[,.(n=.N,histology_pct=100*mean(dc3=='Histology'),strict8140_pct=100*mean(hist==8140),marital_unknown_pct=100*mean(marital3=='Unknown'),race_unknown_pct=100*mean(grepl('Unknown',race,ignore.case=TRUE)|is.na(race)|race=='')),by=year][order(year)]
fwrite(annual,file.path(dest,'50_annual_CUP_composition.csv'))
fwrite(cup[,.(n=.N,percent=100*.N/nrow(cup)),by=stage],file.path(dest,'51_CUP_stage_distribution.csv'))
fwrite(e[keep,.(n=.N,age85plus_n=sum(age>=85),age85plus_pct=100*mean(age>=85),age90plus_n=sum(old90==1)),by=group],file.path(dest,'52_relative_survival_feasibility.csv'))
capture.output(sessionInfo(),file=file.path(dest,'R_sessionInfo.txt'));print(res);print(dg);cat('DONE\n')

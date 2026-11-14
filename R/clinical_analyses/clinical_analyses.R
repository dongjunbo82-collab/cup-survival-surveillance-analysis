options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(splines);library(jsonlite);library(parallel)})
setDTthreads(1)
out<-'outputs/clinical_analyses';dir.create(out,recursive=TRUE,showWarnings=FALSE)
e<-readRDS('work/eligible_cohort.rds');source('R/standardization/analysis_functions.R')
cfg<-fromJSON('config/histology_definitions.json');keep<-e$hist%in%cfg$glandular_codes
metastatic_coverage <- e[keep, .(
  n=.N,
  four_metastatic_sites_complete=sum(organ_complete),
  complete_percent=100*mean(organ_complete)
), by=.(group, era)]
fwrite(metastatic_coverage, file.path(out, 'metastatic_field_audit.csv'))

groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-c('2004-2007','2020-2022')
write_json(list(status='Exploratory clinical analyses',B=2000,seed_rule='20261025+b; shared broad group-era draws',age_groups=c('<65','65-74','>=75'),age_target='Pooled glandular CUP within each age stratum',landmark=60,landmark_endpoint='Survival to diagnosis day 365.25 conditional on survival through day 60; baseline weights, baseline landmark weights',decomposition='Symmetric two-factor decomposition of S365=S90*C; early=(S90late-S90early)*(Clate+Cearly)/2; subsequent=(Clate-Cearly)*(S90late+S90early)/2',inference='Nominal descriptive percentile 95% intervals; descriptive inference'),file.path(out,'specification.json'),pretty=TRUE)
# Use unrounded, refitted and paired existing bootstrap draws for decomposition.
ar<-readRDS('work/glandular_models/summary_arrays.rds')$A
pt<-fread('outputs/glandular_models/15_era_estimates_CI.csv')[definition=='Glandular']
dec<-function(a,b,c,d)c(total=b*d-a*c,early=(b-a)*(d+c)/2,subsequent=(d-c)*(b+a)/2)
dr<-list();dv<-list()
for(g in groups){
 labs<-paste('Era','Glandular',eras,g,sep='__')
 x<-sapply(eras,function(er)c(pt[group==g&period==er&measure=='OS90',estimate],pt[group==g&period==er&measure=='Conditional365given90',estimate]))
 p<-dec(x[1,1],x[1,2],x[2,1],x[2,2]);names(p)<-c('total','early','subsequent');v<-cbind(total=ar[,labs[2],'OS365']-ar[,labs[1],'OS365'],early=(ar[,labs[2],'OS90']-ar[,labs[1],'OS90'])*(ar[,labs[2],'Conditional365given90']+ar[,labs[1],'Conditional365given90'])/2,subsequent=(ar[,labs[2],'Conditional365given90']-ar[,labs[1],'Conditional365given90'])*(ar[,labs[2],'OS90']+ar[,labs[1],'OS90'])/2)
 stopifnot(max(abs(v[,1]-v[,2]-v[,3]))<1e-12,abs(p[1]-p[2]-p[3])<1e-12)
 dv[[g]]<-v
 dr[[g]]<-data.table(group=g,component=names(p),estimate_pp=100*p,lo_pp=100*apply(v,2,quantile,.025),hi_pp=100*apply(v,2,quantile,.975))
}
dd<-rbindlist(dr);for(g in groups[-1]){v<-dv$CUP-dv[[g]];dd<-rbind(dd,data.table(group=paste0('CUP-minus-',g,' gap change'),component=colnames(v),estimate_pp=dd[group=='CUP',estimate_pp]-dd[group==g,estimate_pp],lo_pp=100*apply(v,2,quantile,.025),hi_pp=100*apply(v,2,quantile,.975)))}
fwrite(dd,file.path(out,'decomposition.csv'))
# Refit age-stratum membership models; preserve broad source cohort and shared sampling.
e[,age_stratum:=fcase(age<65,'<65',age<75,'65-74',default='>=75')]
models<-list();points<-list();diagnostics<-list();balances<-list()
for(ag in c('<65','65-74','>=75'))for(er in eras)for(g in groups){
 lab<-paste(ag,er,g,sep='__');tt<-which(keep&e$group=='CUP'&e$age_stratum==ag);ss<-which(keep&e$group==g&as.character(e$era)==er&e$age_stratum==ag)
 m<-build_model(tt,ss,lab,'Age',er,g);f<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(f$error));models[[lab]]<-m
 points[[lab]]<-data.table(label=lab,age=ag,era=er,group=g,measure=measures,estimate=f$values)
 w<-f$odds;diagnostics[[lab]]<-cbind(data.table(label=lab,age=ag,era=er,group=g),as.data.table(as.list(c(f$diagnostic,top1_share=sum(sort(w,decreasing=TRUE)[seq_len(ceiling(length(w)*.01))])/sum(w)))))
 balances[[lab]]<-cbind(data.table(label=lab),f$smd)
}
# Baseline-standardized 60-day sensitivity for all four comparators, first/last eras.
km60<-function(m,w){a<-rowsum(cbind(w,w*m$event),m$timebin,reorder=FALSE);risk<-rev(cumsum(rev(a[,1])));S<-cumprod(pmax(0,1-ifelse(risk>0,a[,2]/risk,0)));j<-which(m$times<=60);s60<-if(length(j))S[max(j)]else 1;j<-which(m$times<=365.25);s365<-if(length(j))S[max(j)]else 1;c(OS60=s60,Conditional365given60=s365/s60)}
base<-list();bp<-list()
for(er in eras)for(g in groups){lab<-paste(er,g,sep='__');m<-build_model(which(keep&e$group=='CUP'),which(keep&e$group==g&as.character(e$era)==er),lab,'Day60',er,g);f<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(f$error));base[[lab]]<-m;bp[[lab]]<-data.table(label=lab,era=er,group=g,measure=c('OS60','Conditional365given60'),estimate=km60(m,f$odds))}
ap<-rbindlist(points);dg<-rbindlist(diagnostics);dg[,support_pass:=n_source>=100&ESS>=200&max_SMD<.1];bpt<-rbindlist(bp)
fwrite(ap,file.path(out,'age_point_estimates.csv'));fwrite(dg,file.path(out,'age_diagnostics.csv'));fwrite(rbindlist(balances),file.path(out,'age_balance.csv'));fwrite(bpt,file.path(out,'day60_point_estimates.csv'))
saveRDS(list(age=models,day60=base),'work/clinical_analyses_models.rds')
repfit<-function(b){set.seed(20261025+b);counts<-integer(N);for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx));v<-matrix(NA_real_,length(models),6);u<-matrix(NA_real_,length(base),2);err<-character()
 for(i in seq_along(models)){f<-tryCatch(fit_model(models[[i]],counts),error=function(ex)list(values=rep(NA_real_,6),error=conditionMessage(ex)));v[i,]<-f$values;if(nzchar(f$error))err<-c(err,paste(names(models)[i],f$error))}
 for(i in seq_along(base)){f<-tryCatch(fit_model(base[[i]],counts,diagnostic=TRUE),error=function(ex)list(error=conditionMessage(ex)));if(nzchar(f$error)){err<-c(err,paste(names(base)[i],f$error))}else u[i,]<-km60(base[[i]],f$odds*counts[base[[i]]$source])}
 list(b=b,v=v,u=u,error=err)}
cp <- 'work/clinical_analyses_bootstrap.rds'
B <- 2000L
signature <- list(
  cohort=unname(tools::md5sum('work/eligible_cohort.rds')),
  code=unname(tools::md5sum('R/clinical_analyses/clinical_analyses.R')),
  functions=unname(tools::md5sum('R/standardization/analysis_functions.R')),
  definitions=unname(tools::md5sum('config/histology_definitions.json')),
  B=B, seed=20261025
)
signature_file <- paste0(cp, '.signature.rds')
if (file.exists(cp)) {
  stopifnot(file.exists(signature_file), identical(signature, readRDS(signature_file)))
}
saveRDS(signature, signature_file)
boot <- if (file.exists(cp)) readRDS(cp) else list()
cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS', '2')),outfile='work/clinical_analyses_workers.log');clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL});clusterExport(cl,c('N','strata','models','base','repfit','fit_model','km_cached','km60'),envir=environment())
for(start in if(length(boot)<B)seq(length(boot)+1L,B,by=50)else integer()){boot<-c(boot,parLapplyLB(cl,start:min(B,start+49L),repfit));saveRDS(boot,cp);cat('Completed',length(boot),'/',B,format(Sys.time()),'\n');flush.console()};stopCluster(cl)
errors<-rbindlist(lapply(boot,function(z)if(length(z$error))data.table(replicate=z$b,error=z$error)else NULL));fwrite(errors,file.path(out,'bootstrap_errors.csv'));stopifnot(nrow(errors)==0)
A<-array(NA_real_,c(B,length(models),6),dimnames=list(NULL,names(models),measures));U<-array(NA_real_,c(B,length(base),2),dimnames=list(NULL,names(base),c('OS60','Conditional365given60')));for(b in 1:B){A[b,,]<-boot[[b]]$v;U[b,,]<-boot[[b]]$u};stopifnot(all(is.finite(A)),all(is.finite(U)))
res<-list();vec<-list()
for(ag in c('<65','65-74','>=75'))for(g in groups){labs<-paste(ag,eras,g,sep='__');for(me in c('OS365','OS90')){sgn<-if(me=='OS90')-1 else 1;p<-sgn*(ap[label==labs[2]&measure==me,estimate]-ap[label==labs[1]&measure==me,estimate]);v<-sgn*(A[,labs[2],me]-A[,labs[1],me]);vec[[paste(ag,g,me)]]<-v;res[[length(res)+1]]<-data.table(age=ag,comparison=g,measure=if(me=='OS90')'Death90_change'else 'OS365_change',estimate_pp=100*p,lo_pp=100*quantile(v,.025),hi_pp=100*quantile(v,.975))}}
rr<-rbindlist(res);gaps<-list()
for(ag in c('<65','65-74','>=75'))for(g in groups[-1]){v<-vec[[paste(ag,'CUP','OS365')]]-vec[[paste(ag,g,'OS365')]];p<-rr[age==ag&comparison=='CUP'&measure=='OS365_change',estimate_pp]-rr[age==ag&comparison==g&measure=='OS365_change',estimate_pp];vec[[paste(ag,g,'gap')]]<-v;gaps[[length(gaps)+1]]<-data.table(age=ag,control=g,estimate_pp=p,lo_pp=100*quantile(v,.025),hi_pp=100*quantile(v,.975),support_pass=all(dg[age==ag&group%in%c('CUP',g),support_pass]))}
gg<-rbindlist(gaps);fwrite(rr,file.path(out,'age_changes.csv'));fwrite(gg,file.path(out,'age_gap_changes.csv'))
het<-list();pairs<-combn(c('<65','65-74','>=75'),2,simplify=FALSE)
for(pa in pairs)for(g in groups[-1]){v<-vec[[paste(pa[2],g,'gap')]]-vec[[paste(pa[1],g,'gap')]];het[[length(het)+1]]<-data.table(age_contrast=paste(pa[2],'-',pa[1]),control=g,estimate_pp=gg[age==pa[2]&control==g,estimate_pp]-gg[age==pa[1]&control==g,estimate_pp],lo_pp=100*quantile(v,.025),hi_pp=100*quantile(v,.975))}
fwrite(rbindlist(het),file.path(out,'age_between_strata.csv'))
ur<-list();for(er in eras)for(me in c('OS60','Conditional365given60'))for(g in groups){lab<-paste(er,g,sep='__');v<-U[,lab,me];p<-bpt[label==lab&measure==me,estimate];ur[[length(ur)+1]]<-data.table(era=er,comparison=g,measure=me,estimate_pp=100*p,lo_pp=100*quantile(v,.025),hi_pp=100*quantile(v,.975));if(g!='CUP'){v<-U[,paste(er,'CUP',sep='__'),me]-v;p<-bpt[label==paste(er,'CUP',sep='__')&measure==me,estimate]-p;ur[[length(ur)+1]]<-data.table(era=er,comparison=paste('CUP-minus',g),measure=me,estimate_pp=100*p,lo_pp=100*quantile(v,.025),hi_pp=100*quantile(v,.975))}}
fwrite(rbindlist(ur),file.path(out,'day60_estimates_CI.csv'));saveRDS(list(A=A,U=U),'work/clinical_analyses_arrays.rds')
write_json(list(B=B,all_finite=TRUE,failed_fits=0,age_models=nrow(dg),age_models_supported=sum(dg$support_pass),max_age_SMD=max(dg$max_SMD),min_age_ESS=min(dg$ESS)),file.path(out,'run_QA.json'),pretty=TRUE);capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'));cat('DONE\n')

# Independent weighted survival checks across age and boundary models.
survival_checks <- list()
for (kind in c('age', 'day60')) {
  selected_models <- if (kind == 'age') models else base
  for (label in names(selected_models)) {
    model <- selected_models[[label]]
    fit <- fit_model(model, diagnostic=TRUE)
    curve <- survival::survfit(
      survival::Surv(e$t[model$source], e$event[model$source]) ~ 1,
      weights=fit$odds
    )
    independent <- summary(curve, times=c(60, 90, 365.25), extend=TRUE)$surv
    cached <- km_cached(model, fit$odds)
    difference <- max(abs(cached[c('OS90', 'OS365')] - independent[2:3]))
    if (kind == 'day60') {
      boundary <- km60(model, fit$odds)
      difference <- max(difference, abs(boundary['OS60'] - independent[1]))
    }
    stopifnot(difference < 1e-10)
    survival_checks[[length(survival_checks)+1L]] <- data.table(
      kind=kind, label=label, max_absolute_difference=difference
    )
  }
}
fwrite(rbindlist(survival_checks), file.path(out, 'independent_survfit_QA.csv'))

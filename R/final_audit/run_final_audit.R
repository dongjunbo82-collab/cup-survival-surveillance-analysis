# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
options(stringsAsFactors=FALSE)
Sys.setenv(STAGE3_FUNCTIONS_ONLY='1')
source('R/validation_stage3/run_validation.R',encoding='UTF-8')
Sys.unsetenv('STAGE3_FUNCTIONS_ONLY')
out<-file.path(root,'outputs/final_audit');dir.create(out,recursive=TRUE,showWarnings=FALSE)
groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-levels(e$era)
write_tab<-function(name,z)fwrite(z,file.path(out,paste0(name,'.csv')),bom=TRUE,na='NA')
models<-list();points<-list();diags<-list();balance<-list();curves<-list();risks<-list();support<-list()
full_curve<-function(m,w){
 a<-rowsum(cbind(w,w*m$event),m$timebin,reorder=FALSE)
 risk<-rev(cumsum(rev(a[,1])));S<-cumprod(pmax(0,1-a[,2]/risk))
 data.table(day=c(0,m$times[m$times<=365.25],365.25),survival=c(1,S[m$times<=365.25],if(any(m$times<=365.25))tail(S[m$times<=365.25],1)else 1))
}
add_model<-function(def,period,g,annual=FALSE){
 target<-which(e$group=='CUP'&(def=='Broad'|e$hist==8140))
 source<-which(e$group==g&(if(annual)e$year==as.numeric(period)else as.character(e$era)==period)&(def=='Broad'|e$hist==8140))
 lab<-paste(if(annual)'Annual'else'Era',def,period,g,sep='__')
 m<-build_model(target,source,lab,def,period,g);f<-fit_model(m,diagnostic=TRUE)
 if(nzchar(f$error))stop(lab,':',f$error)
 models[[lab]]<<-m
 points[[lab]]<<-data.table(label=lab,scope=if(annual)'Annual'else'Era',definition=def,period=as.character(period),group=g,measure=measures,estimate=f$values)
 w<-f$odds;tw<-sum(w);ord<-order(w,decreasing=TRUE);top<-head(ord,max(1,ceiling(length(w)*.01)))
 d<-c(f$diagnostic,weight_p50=unname(quantile(w,.5)),weight_p95=unname(quantile(w,.95)),weight_p99=unname(quantile(w,.99)),top1pct_weight_share=sum(w[top])/tw)
 diags[[lab]]<<-cbind(data.table(label=lab,scope=if(annual)'Annual'else'Era',definition=def,period=as.character(period),group=g),as.data.table(as.list(d)))
 balance[[lab]]<<-cbind(data.table(label=lab),f$smd)
 # Empirical covariate support: source range and absence of each target category.
 for(v in covnames){
  tv<-e[[v]][target];sv<-e[[v]][source]
  outside<-if(v=='age')mean(tv<min(sv)|tv>max(sv))else mean(!tv%in%unique(sv))
  support[[paste(lab,v)]]<<-data.table(label=lab,variable=v,target_outside_source_support_pct=100*outside)
 }
 if(!annual){
  curves[[lab]]<<-cbind(data.table(label=lab,definition=def,era=period,group=g),full_curve(m,w))
  risks[[lab]]<<-rbindlist(lapply(c(0,30,90,180,365.25),function(tau){
   keep<-e$t[m$source]>=tau;ww<-w[keep]
   data.table(label=lab,definition=def,era=period,group=g,day=tau,raw_at_risk=sum(keep),weighted_at_risk_normalized=sum(ww)/tw*length(w),risk_ESS=if(length(ww))sum(ww)^2/sum(ww^2)else 0)
  }))
 }
}
for(def in c('Broad','Strict8140')){
 for(er in eras)for(g in groups)add_model(def,er,g)
 for(yr in 2004:2022)for(g in groups)add_model(def,as.character(yr),g,TRUE)
}
pt<-rbindlist(points);dg<-rbindlist(diags)
dg[,support_pass:=n_source>=100&ESS>=200&max_SMD<.10]
write_tab('01_point_estimates',pt);write_tab('02_model_diagnostics',dg)
write_tab('03_model_balance',rbindlist(balance));write_tab('04_covariate_support',rbindlist(support))
write_tab('05_weighted_curves',rbindlist(curves));write_tab('06_curve_risk_support',rbindlist(risks))
write_tab('07_yearly_code_drift',e[,.(n=.N,strict8140_pct=100*mean(hist==8140),histology_pct=100*mean(dc3=='Histology'),cytology_pct=100*mean(dc3=='Cytology'),micro_unspecified_pct=100*mean(dc3=='Micro_unspecified'),complete_dates_pct=100*mean(grepl('^Complete',flag)),alive_censor_under1y_pct=100*mean(event==0&t<365.25),organ_complete_pct=100*mean(organ_complete)),by=.(group,year)][order(group,year)])
write_tab('08_era_histology_composition',e[,.N,by=.(group,era,hist)][,percent:=100*N/sum(N),by=.(group,era)][order(group,era,-N)])
write_tab('09_era_covariate_composition',rbindlist(lapply(c('sex','race','marital3','dc3','age_band','income4','rural3'),function(v)e[,.N,by=c('group','era',v)][,.(group,era,variable=v,category=get(v),n=N)][,percent:=100*n/sum(n),by=.(group,era)])))
write_tab('10_era_followup',e[,.(n=.N,dead=sum(event),zero_days=sum(days==0),alive_zero=sum(days==0&event==0),zero_months=sum(months==0),alive_censor_under1y=sum(event==0&t<365.25),incomplete_dates=sum(!grepl('^Complete',flag))),by=.(group,era)])
write_tab('11_yearly_organ_labels',melt(e[,.(group,year,bone,brain,liver,lung)],id.vars=c('group','year'),variable.name='organ',value.name='label')[,.N,by=.(group,year,organ,label)])
write_tab('12_year_composition',e[,.N,by=.(group,era,year)][,proportion:=N/sum(N),by=.(group,era)])
saveRDS(models,file.path(root,'work/final_audit_models.rds'))
cat('Models:',length(models),' Annual support failures:',nrow(dg[scope=='Annual'&!support_pass]),'\n')
print(dg[scope=='Annual'&!support_pass,.(definition,period,group,n_source,ESS,max_SMD)])
flush.console()
# Pair the additional models with the original seed and patient draws; original
# 30-era-model replicates can then be reused without resampling inconsistency.
bootnames<-names(models)[vapply(models,function(m)!grepl('^Era__',m$label)||m$era=='2008-2014',logical(1))]
newmodels<-models[bootnames]
replicate_final<-function(b){
 set.seed(20261025+b);counts<-integer(N)
 for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx))
 v<-matrix(NA_real_,length(newmodels),length(measures));err<-character(length(newmodels))
 for(i in seq_along(newmodels)){
  f<-tryCatch(fit_model(newmodels[[i]],counts),error=function(ex)list(values=rep(NA_real_,6),error=conditionMessage(ex)))
  v[i,]<-f$values;err[i]<-f$error
 }
 # q resampled within original group-era strata, sharing reference CUP draws.
 q<-as.numeric(rowsum(counts[cupidx],cupyears,reorder=TRUE))
 qs<-as.numeric(rowsum(counts[strictidx],strictyears,reorder=TRUE))
 list(b=b,v=v,error=err,q=q,qs=qs)
}
e_year<-e$year;e_group<-e$group;e_hist<-e$hist
cupidx<-which(e_group=='CUP');strictidx<-which(e_group=='CUP'&e_hist==8140)
cupyears<-e_year[cupidx];strictyears<-e_year[strictidx]
if('--points-only'%in%commandArgs(TRUE))quit(save='no')
if('--bench'%in%commandArgs(TRUE)){
 print(system.time(z<-lapply(1:3,replicate_final)));print(lapply(z,function(z)table(z$error)));quit(save='no')
}
B<-2000L
signature<-list(input_md5=unname(tools::md5sum(file.path(root,'work/stage2_eligible.rds'))),labels=bootnames,seed=20261025,B=B)
sigfile<-file.path(root,'work/final_audit_signature.rds');checkpoint<-file.path(root,'work/final_audit_checkpoint.rds')
if(file.exists(sigfile))stopifnot(identical(signature,readRDS(sigfile)))else saveRDS(signature,sigfile)
allboot<-if(file.exists(checkpoint))readRDS(checkpoint)else list()
stopifnot(length(allboot)<=B)
cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS','2')),outfile=file.path(root,'work/final_audit_workers.log'))
clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL})
clusterExport(cl,c('N','strata','newmodels','measures','replicate_final','fit_model','km_cached','cupidx','strictidx','cupyears','strictyears'),envir=environment())
starts<-if(length(allboot)<B)seq(length(allboot)+1L,B,by=25)else integer()
for(start in starts){
 z<-parLapplyLB(cl,start:min(B,start+24),replicate_final)
 allboot<-c(allboot,z);saveRDS(allboot,checkpoint)
 cat('Completed',length(allboot),'/',B,'at',format(Sys.time()),'\n');flush.console()
}
stopCluster(cl)
saveRDS(allboot,file.path(out,'bootstrap_replicates.rds'))
write_tab('13_bootstrap_failures',rbindlist(lapply(allboot,function(z)data.table(replicate=z$b,label=bootnames,error=z$error)))[nzchar(error)])
saveRDS(bootnames,file.path(out,'bootstrap_labels.rds'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
writeLines('Additional bootstrap finished; summaries and independent numerical checks required.',file.path(out,'BOOTSTRAP_COMPLETE.txt'))

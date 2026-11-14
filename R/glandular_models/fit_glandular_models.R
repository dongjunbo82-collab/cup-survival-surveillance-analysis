# Generated artifacts are written to local output directories.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/glandular_models', 'work/landmark',
            'outputs/standardization', 'outputs/cohort_checks',
            'outputs/definition_models', 'outputs/era_models', 'outputs/glandular_models',
            'outputs/pandemic_sensitivity', 'outputs/absolute_outcomes',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
Sys.setenv(CUP_FUNCTIONS_ONLY='1')
source('R/definition_models/fit_definition_models.R',encoding='UTF-8')
Sys.unsetenv('CUP_FUNCTIONS_ONLY');library(jsonlite);setDTthreads(1)
out<-file.path(root,'outputs/glandular_models');cfg<-fromJSON(file.path(root,'config/histology_definitions.json'))
groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-levels(e$era)
wt<-function(n,x)fwrite(x,file.path(out,paste0(n,'.csv')),bom=TRUE,na='NA')
keep<-list(Glandular=e$hist%in%cfg$glandular_codes,MinimalNoNE=e$hist%in%cfg$minimal_no_NE_keep_codes)
models<-list();pts<-list();dgs<-list();balances<-list();curves<-list();risks<-list();references<-list()
for(def in names(keep)){
 target<-which(keep[[def]]&e$group=='CUP')
 for(scope in if(def=='Glandular')c('Era','Annual')else'Era')for(period in if(scope=='Era')eras else as.character(2004:2022))for(g in groups){
  source<-which(keep[[def]]&e$group==g&(if(scope=='Era')as.character(e$era)==period else e$year==as.integer(period)))
  label<-paste(scope,def,period,g,sep='__');m<-build_model(target,source,label,def,period,g);f<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(f$error))
  if(is.null(references[[def]]))references[[def]]<-list(target=m$target,bin=m$targetbin)
  stopifnot(identical(m$targetbin,references[[def]]$bin))
  models[[label]]<-m;pts[[label]]<-data.table(label=label,scope=scope,definition=def,period=period,group=g,measure=measures,estimate=f$values)
  w<-f$odds;dgs[[label]]<-cbind(data.table(label=label,scope=scope,definition=def,period=period,group=g),as.data.table(as.list(c(f$diagnostic,weight_p99=unname(quantile(w,.99)),top1pct_share=sum(sort(w,decreasing=TRUE)[seq_len(max(1,ceiling(length(w)*.01)))])/sum(w)))))
  balances[[label]]<-cbind(data.table(label=label),f$smd)
  if(scope=='Era'&def=='Glandular'){
   a<-rowsum(cbind(w,w*m$event),m$timebin,reorder=FALSE);s<-cumprod(pmax(0,1-a[,2]/rev(cumsum(rev(a[,1])))))
   curves[[label]]<-data.table(era=period,group=g,day=c(0,m$times[m$times<=365.25],365.25),survival=c(1,s[m$times<=365.25],tail(s[m$times<=365.25],1)))
   risks[[label]]<-rbindlist(lapply(c(0,30,90,180,365.25),function(tau){jj<-e$t[m$source]>=tau;ww<-w[jj];data.table(era=period,group=g,day=tau,n_risk=sum(jj),weighted_risk=sum(ww)/sum(w)*length(w),risk_ESS=sum(ww)^2/sum(ww^2))}))
  }
 }
}
pt<-rbindlist(pts);dg<-rbindlist(dgs);dg[,support_pass:=n_source>=100&ESS>=200&max_SMD<.1]
wt('04_point_estimates',pt);wt('05_weight_diagnostics',dg);wt('06_balance',rbindlist(balances));wt('07_glandular_curves',rbindlist(curves));wt('08_risk_support',rbindlist(risks))
wt('09_glandular_year_composition',e[keep$Glandular,.N,by=.(group,era,year)][,proportion:=N/sum(N),by=.(group,era)])
wt('10_glandular_baseline',rbindlist(lapply(c('sex','race','age_band','marital3','dc3'),function(v)e[keep$Glandular,.N,by=c('group',v)][,.(group,variable=v,category=get(v),n=N)][,percent:=100*n/sum(n),by=group])))
wt('11_glandular_histology_counts',e[keep$Glandular,.N,by=.(group,era,hist)])
wt('12_excluded_by_code_year',e[!keep$Glandular,.N,by=.(group,year,hist)])
saveRDS(models,file.path(root,'work/glandular_models/models.rds'))
fit_fast<-function(m,counts,refcounts){
 nt<-numeric(m$nb);nt[seq_along(refcounts)]<-refcounts
 rs<-rowsum(counts[m$source],m$sourcebin,reorder=TRUE);ns<-numeric(m$nb);ns[as.integer(rownames(rs))]<-rs[,1]
 total<-nt+ns;y<-ifelse(total>0,nt/pmax(total,1),.5)
 f<-suppressWarnings(glm.fit(m$X,y,weights=total,family=binomial(),control=glm.control(maxit=40)))
 if(!f$converged)return(list(values=rep(NA_real_,6),error='nonconvergence'))
 p<-f$fitted.values[m$sourcebin]
 if(any(!is.finite(p)|p>=1|p<=0))return(list(values=rep(NA_real_,6),error='boundary_probability'))
 list(values=km_cached(m,p/(1-p)*counts[m$source]),error='')
}
cupidx<-references$Glandular$target;cupyear<-e$year[cupidx]
repfit<-function(b){
 set.seed(20261025+b);counts<-integer(N)
 for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx))
 refs<-lapply(references,function(x)as.numeric(rowsum(counts[x$target],x$bin,reorder=TRUE)))
 v<-matrix(NA_real_,length(models),6);err<-character(length(models))
 for(i in seq_along(models)){m<-models[[i]];z<-tryCatch(fit_fast(m,counts,refs[[m$definition]]),error=function(ex)list(values=rep(NA_real_,6),error=conditionMessage(ex)));v[i,]<-z$values;err[i]<-z$error}
 list(b=b,v=v,error=err,q=as.numeric(rowsum(counts[cupidx],cupyear,reorder=TRUE)))
}
# Optimization only caches repeated reference counts; verify against original code.
set.seed(20261026);counts<-integer(N)
for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx))
qa<-list()
for(lab in c('Era__Glandular__2020-2022__Pancreas','Annual__Glandular__2014__Stomach','Era__MinimalNoNE__2004-2007__CUP')){
 m<-models[[lab]];ref<-as.numeric(rowsum(counts[m$target],m$targetbin,reorder=TRUE));a<-fit_fast(m,counts,ref);b<-fit_model(m,counts)
 error<-max(abs(a$values-b$values));stopifnot(error<1e-10);qa[[lab]]<-data.table(check='Cached vs original reference likelihood',label=lab,error=error,passed=TRUE)
}
wt('13_optimization_QA',rbindlist(qa));cat('Point models:',length(models),'Annual failures:',nrow(dg[scope=='Annual'&support_pass==FALSE]),'\n');flush.console()
if('--bench'%in%commandArgs(TRUE)){print(system.time(z<-lapply(1:3,repfit)));quit(save='no')}
B<-2000L;cp<-file.path(root,'work/glandular_models/bootstrap_checkpoint.rds')
sig<-list(input=unname(tools::md5sum(file.path(root,'work/eligible_cohort.rds'))),cfg=unname(tools::md5sum(file.path(root,'config/histology_definitions.json'))),labels=names(models),B=B)
sf<-file.path(root,'work/glandular_models/signature.rds');if(file.exists(sf))stopifnot(identical(sig,readRDS(sf)))else saveRDS(sig,sf)
boot<-if(file.exists(cp))readRDS(cp)else list()
cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS','2')),outfile=file.path(root,'work/glandular_models/workers.log'))
clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL})
clusterExport(cl,c('N','strata','references','models','repfit','fit_fast','km_cached','cupidx','cupyear'),envir=environment())
for(start in if(length(boot)<B)seq(length(boot)+1L,B,by=60)else integer()){
 boot<-c(boot,parLapplyLB(cl,start:min(B,start+59L),repfit));saveRDS(boot,cp)
 cat('Completed',length(boot),'/',B,format(Sys.time()),'\n');flush.console()
}
stopCluster(cl);saveRDS(boot,file.path(out,'bootstrap_replicates.rds'));saveRDS(names(models),file.path(out,'bootstrap_labels.rds'))
wt('14_bootstrap_failures',rbindlist(lapply(boot,function(z)data.table(replicate=z$b,label=names(models),error=z$error)))[nzchar(error)])
stopifnot(mean(vapply(boot,function(z)all(z$error==''),logical(1)))>=.95)
writeLines('2000 completed; summaries required',file.path(out,'BOOTSTRAP_COMPLETE.txt'));capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))

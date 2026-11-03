# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
Sys.setenv(STAGE3_FUNCTIONS_ONLY='1')
source('R/validation_stage3/run_validation.R',encoding='UTF-8')
Sys.unsetenv('STAGE3_FUNCTIONS_ONLY');library(jsonlite);setDTthreads(1)
out<-file.path(root,'outputs/histology_revision');cfg<-fromJSON(file.path(root,'config/histology_definitions.json'))
groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-levels(e$era)
wt<-function(n,x)fwrite(x,file.path(out,paste0(n,'.csv')),bom=TRUE,na='NA')
keep<-list(Reviewed=e$hist%in%cfg$reviewed_keep_codes,MinimalNoNE=e$hist%in%cfg$minimal_no_NE_keep_codes)
models<-list();pts<-list();dgs<-list();balances<-list();curves<-list();risks<-list();references<-list()
for(def in names(keep)){
 target<-which(keep[[def]]&e$group=='CUP')
 for(scope in if(def=='Reviewed')c('Era','Annual')else'Era')for(period in if(scope=='Era')eras else as.character(2004:2022))for(g in groups){
  source<-which(keep[[def]]&e$group==g&(if(scope=='Era')as.character(e$era)==period else e$year==as.integer(period)))
  label<-paste(scope,def,period,g,sep='__');m<-build_model(target,source,label,def,period,g);f<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(f$error))
  if(is.null(references[[def]]))references[[def]]<-list(target=m$target,bin=m$targetbin)
  stopifnot(identical(m$targetbin,references[[def]]$bin))
  models[[label]]<-m;pts[[label]]<-data.table(label=label,scope=scope,definition=def,period=period,group=g,measure=measures,estimate=f$values)
  w<-f$odds;dgs[[label]]<-cbind(data.table(label=label,scope=scope,definition=def,period=period,group=g),as.data.table(as.list(c(f$diagnostic,weight_p99=unname(quantile(w,.99)),top1pct_share=sum(sort(w,decreasing=TRUE)[seq_len(max(1,ceiling(length(w)*.01)))])/sum(w)))))
  balances[[label]]<-cbind(data.table(label=label),f$smd)
  if(scope=='Era'&def=='Reviewed'){
   a<-rowsum(cbind(w,w*m$event),m$timebin,reorder=FALSE);s<-cumprod(pmax(0,1-a[,2]/rev(cumsum(rev(a[,1])))))
   curves[[label]]<-data.table(era=period,group=g,day=c(0,m$times[m$times<=365.25],365.25),survival=c(1,s[m$times<=365.25],tail(s[m$times<=365.25],1)))
   risks[[label]]<-rbindlist(lapply(c(0,30,90,180,365.25),function(tau){jj<-e$t[m$source]>=tau;ww<-w[jj];data.table(era=period,group=g,day=tau,n_risk=sum(jj),weighted_risk=sum(ww)/sum(w)*length(w),risk_ESS=sum(ww)^2/sum(ww^2))}))
  }
 }
}
pt<-rbindlist(pts);dg<-rbindlist(dgs);dg[,support_pass:=n_source>=100&ESS>=200&max_SMD<.1]
wt('04_point_estimates',pt);wt('05_weight_diagnostics',dg);wt('06_balance',rbindlist(balances));wt('07_reviewed_curves',rbindlist(curves));wt('08_risk_support',rbindlist(risks))
wt('09_reviewed_year_composition',e[keep$Reviewed,.N,by=.(group,era,year)][,proportion:=N/sum(N),by=.(group,era)])
wt('10_reviewed_baseline',rbindlist(lapply(c('sex','race','age_band','marital3','dc3'),function(v)e[keep$Reviewed,.N,by=c('group',v)][,.(group,variable=v,category=get(v),n=N)][,percent:=100*n/sum(n),by=group])))
wt('11_reviewed_histology_counts',e[keep$Reviewed,.N,by=.(group,era,hist)])
wt('12_excluded_by_code_year',e[!keep$Reviewed,.N,by=.(group,year,hist)])
saveRDS(models,file.path(root,'work/histology_review/models.rds'))
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
cupidx<-references$Reviewed$target;cupyear<-e$year[cupidx]
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
for(lab in c('Era__Reviewed__2020-2022__Pancreas','Annual__Reviewed__2014__Stomach','Era__MinimalNoNE__2004-2007__CUP')){
 m<-models[[lab]];ref<-as.numeric(rowsum(counts[m$target],m$targetbin,reorder=TRUE));a<-fit_fast(m,counts,ref);b<-fit_model(m,counts)
 error<-max(abs(a$values-b$values));stopifnot(error<1e-10);qa[[lab]]<-data.table(check='Cached vs original reference likelihood',label=lab,error=error,passed=TRUE)
}
wt('13_optimization_QA',rbindlist(qa));cat('Point models:',length(models),'Annual failures:',nrow(dg[scope=='Annual'&support_pass==FALSE]),'\n');flush.console()
if('--bench'%in%commandArgs(TRUE)){print(system.time(z<-lapply(1:3,repfit)));quit(save='no')}
B<-2000L;cp<-file.path(root,'work/histology_review/bootstrap_checkpoint.rds')
sig<-list(input=unname(tools::md5sum(file.path(root,'work/stage2_eligible.rds'))),cfg=unname(tools::md5sum(file.path(root,'config/histology_definitions.json'))),labels=names(models),B=B)
sf<-file.path(root,'work/histology_review/signature.rds');if(file.exists(sf))stopifnot(identical(sig,readRDS(sf)))else saveRDS(sig,sf)
boot<-if(file.exists(cp))readRDS(cp)else list()
cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS','2')),outfile=file.path(root,'work/histology_review/workers.log'))
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

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
suppressPackageStartupMessages({library(data.table);library(splines);library(parallel)})
setDTthreads(2)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'))
out<-file.path(root,'outputs/definition_models');dir.create(out,recursive=TRUE,showWarnings=FALSE)
e<-readRDS(file.path(root,'work/eligible_cohort.rds'))
stopifnot(nrow(e)==381979,uniqueN(e$id)==nrow(e))
covnames<-c('age','old90','sex','race','marital3','dc3')
e[,key:=.GRP,by=covnames]
N<-nrow(e)
strata<-split(seq_len(N),interaction(e$group,e$era,drop=TRUE))
measures<-c('OS90','OS365','RMST90','RMST365','RMST90to365','Conditional365given90')
km_cached<-function(m,w){
 a<-rowsum(cbind(w,w*m$event),m$timebin,reorder=FALSE)
 risk<-rev(cumsum(rev(a[,1])))
 hz<-ifelse(risk>0,a[,2]/risk,0)
 S<-cumprod(pmax(0,1-hz));tm<-m$times
 sval<-function(tau){j<-which(tm<=tau);if(length(j))S[max(j)]else 1}
 area<-function(tau){j<-which(tm<tau);sum(diff(c(0,tm[j],tau))*c(1,S[j]))}
 s90<-sval(90);s365<-sval(365.25);a90<-area(90);a365<-area(365.25)
 c(OS90=s90,OS365=s365,RMST90=a90,RMST365=a365,RMST90to365=a365-a90,Conditional365given90=if(s90>0)s365/s90 else NA_real_)
}
build_model<-function(target,source,label,definition,era,g,extra=NULL){
 source<-source[order(e$t[source])]
 zz<-rbind(e[target,covnames,with=FALSE],e[source,covnames,with=FALSE])
 tr<-c('ns(age,df=4)','old90','factor(sex)','factor(race)','factor(marital3)','factor(dc3)')
 tr<-tr[vapply(tr,function(v){if(grepl('^factor',v))uniqueN(zz[[sub('^factor\\((.*)\\)$','\\1',v)]])>1 else if(v=='old90')uniqueN(zz$old90)>1 else TRUE},logical(1))]
 Xfull<-model.matrix(as.formula(paste('~',paste(tr,collapse='+'))),data=zz)
 keys<-c(e$key[target],e$key[source]);uk<-unique(keys);bins<-match(keys,uk)
 X<-Xfull[match(uk,keys),,drop=FALSE]
 nt<-tabulate(head(bins,length(target)),nbins=length(uk));ns<-tabulate(tail(bins,length(source)),nbins=length(uk))
 times<-sort(unique(e$t[source]));timebin<-match(e$t[source],times)
 list(label=label,definition=definition,era=era,group=g,target=target,source=source,X=X,nb=length(uk),targetbin=head(bins,length(target)),sourcebin=tail(bins,length(source)),nt=nt,ns=ns,event=e$event[source],times=times,timebin=timebin,source_days=e$days[source],source_months=e$months[source])
}
fit_model<-function(m,counts=NULL,diagnostic=FALSE){
 if(is.null(counts)){nt<-m$nt;ns<-m$ns;sc<-rep(1,length(m$source))}else{
  # Reference may not cover every union key; fill via indexed rowsum labels.
  rt<-rowsum(counts[m$target],m$targetbin,reorder=TRUE);nt<-numeric(m$nb);nt[as.integer(rownames(rt))]<-rt[,1]
  rs<-rowsum(counts[m$source],m$sourcebin,reorder=TRUE);ns<-numeric(m$nb);ns[as.integer(rownames(rs))]<-rs[,1]
  sc<-counts[m$source]
 }
 total<-nt+ns;y<-ifelse(total>0,nt/pmax(total,1),0.5)
 f<-suppressWarnings(glm.fit(x=m$X,y=y,weights=total,family=binomial(),control=glm.control(maxit=40)))
 if(!f$converged)return(list(values=rep(NA_real_,6),error='nonconvergence'))
 p<-f$fitted.values[m$sourcebin]
 if(any(!is.finite(p)|p>=1|p<=0))return(list(values=rep(NA_real_,6),error='boundary_probability'))
 odds<-p/(1-p);w<-odds*sc
 v<-km_cached(m,w)
 ans<-list(values=v,error='')
 if(diagnostic){
  X<-m$X[,-1,drop=FALSE]
  mt<-colSums(X*nt)/sum(nt);ms<-colSums(X*ns)/sum(ns)
  vt<-(colSums(X^2*nt)/sum(nt)-mt^2)*sum(nt)/(sum(nt)-1)
  vs<-(colSums(X^2*ns)/sum(ns)-ms^2)*sum(ns)/(sum(ns)-1)
  wb<-numeric(m$nb);rs<-rowsum(w,m$sourcebin,reorder=TRUE);wb[as.integer(rownames(rs))]<-rs[,1]
  mw<-colSums(X*wb)/sum(wb);den<-sqrt(pmax(0,(vt+vs)/2));smd<-(mt-mw)/den
  ans$diagnostic<-c(n_target=sum(nt),n_source=sum(ns),ESS=sum(w)^2/sum(odds^2*sc),weight_max=max(odds[sc>0]),max_SMD=max(abs(smd),na.rm=TRUE))
  ans$smd<-data.table(variable=colnames(X),SMD=smd)
  ans$odds<-odds
 }
 ans
}
if(Sys.getenv('CUP_FUNCTIONS_ONLY')!='1'){
models<-list();point<-list();diagnostics<-list();balance<-list()
for(def in c('Broad','Strict8140')){
 target<-which(e$group=='CUP'&(def=='Broad'|e$hist==8140))
 for(er in c('2004-2007','2015-2019','2020-2022'))for(g in c('CUP','Lung','Pancreas','Colorectal','Stomach')){
  source<-which(e$group==g&as.character(e$era)==er&(def=='Broad'|e$hist==8140))
  label<-paste(def,er,g,sep='__');m<-build_model(target,source,label,def,er,g);f<-fit_model(m,diagnostic=TRUE)
  if(nzchar(f$error))stop(label,':',f$error)
  models[[label]]<-m
  point[[label]]<-data.table(label=label,definition=def,era=er,group=g,measure=measures,estimate=f$values)
  diagnostics[[label]]<-cbind(data.table(label=label,definition=def,era=er,group=g),as.data.table(as.list(f$diagnostic)))
  balance[[label]]<-cbind(data.table(label=label),f$smd)
 }
}
pt<-rbindlist(point);dg<-rbindlist(diagnostics)
fwrite(pt,file.path(out,'01_fixed_target_point_estimates.csv'),bom=TRUE)
fwrite(dg,file.path(out,'02_weight_diagnostics.csv'),bom=TRUE)
fwrite(rbindlist(balance),file.path(out,'03_balance.csv'),bom=TRUE)
# Compressed binomial likelihood must reproduce the full-record glm estimates.
reference_result<-fread(file.path(root,'outputs/cohort_checks/17_fixed_target_era_OS.csv'))
chk<-merge(pt[measure=='OS365'],reference_result[,.(definition,era,group,reference_result=OS365)],by=c('definition','era','group'))
stopifnot(max(abs(chk$estimate-chk$reference_result))<1e-7)
saveRDS(models,file.path(root,'work/definition_model_cache.rds'))
replicate_fit<-function(b){
 set.seed(20261025+b)
 counts<-integer(N)
 for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx))
 v<-matrix(NA_real_,length(models),length(measures));err<-character(length(models))
 for(i in seq_along(models)){
  ans<-tryCatch(fit_model(models[[i]],counts),error=function(ex)list(values=rep(NA_real_,6),error=conditionMessage(ex)))
  v[i,]<-ans$values;err[i]<-ans$error
 }
 list(b=b,v=v,error=err)
}
if('--bench'%in%commandArgs(trailingOnly=TRUE)){
 print(system.time(z<-lapply(1:3,replicate_fit)))
 print(lapply(z,function(x)table(x$error)))
 quit(save='no')
}
B<-2000L
cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS','2')),outfile=file.path(root,'work/definition_workers.log'))
clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL})
clusterExport(cl,c('N','strata','models','measures','replicate_fit','fit_model','km_cached'),envir=environment())
checkpoint<-file.path(root,'work/definition_bootstrap_checkpoint.rds')
signature<-list(input_md5=unname(tools::md5sum(file.path(root,'work/eligible_cohort.rds'))),labels=names(models),measures=measures,seed_base=20261025,N=N)
signature_file<-file.path(root,'work/definition_checkpoint_signature.rds')
if(file.exists(signature_file))stopifnot(identical(signature,readRDS(signature_file)))else saveRDS(signature,signature_file)
allboot<-if(file.exists(checkpoint))readRDS(checkpoint)else list()
stopifnot(length(allboot)<=B)
if(length(allboot))stopifnot(identical(vapply(allboot,function(x)as.integer(x$b),integer(1)),seq_len(length(allboot))))
starts<-if(length(allboot)<B)seq(length(allboot)+1L,B,by=50)else integer(0)
for(start in starts){
 ids<-start:min(B,start+49)
 x<-parLapplyLB(cl,ids,replicate_fit)
 allboot<-c(allboot,x)
 saveRDS(allboot,file.path(root,'work/definition_bootstrap_checkpoint.rds'))
 cat('Completed',length(allboot),'/',B,'replicates at',format(Sys.time()),'\n');flush.console()
}
stopCluster(cl)
saveRDS(allboot,file.path(out,'bootstrap_replicates.rds'))
failure<-rbindlist(lapply(allboot,function(x)data.table(replicate=x$b,label=names(models),error=x$error)))
fwrite(failure[nzchar(error)],file.path(out,'04_bootstrap_failures.csv'),bom=TRUE)
if(mean(vapply(allboot,function(x)all(!nzchar(x$error)),logical(1)))<.95)stop('Fewer than95% all-model successful replicates; inspect failures')
writeLines('Bootstrap completed; next run summarize_definition_models.R.',file.path(out,'BOOTSTRAP_COMPLETE.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
}

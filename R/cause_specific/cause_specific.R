library(data.table)
library(splines)
library(parallel)
library(jsonlite)
setDTthreads(1)
out <- 'outputs/cause_specific'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
e <- readRDS('work/eligible_cohort.rds')
source('R/standardization/analysis_functions.R')
cfg <- fromJSON('config/histology_definitions.json')
keep <- e$hist %in% cfg$glandular_codes
groups <- c('CUP','Lung','Pancreas','Colorectal','Stomach')
eras <- c('2004-2007','2020-2022')
src <- Sys.getenv('CUP_CSS_EXPORT_DIR')
stopifnot(nzchar(src), dir.exists(src))
stopifnot(sum(keep)==380648)
e[,css:=NA_integer_]
for(g in groups){
 f <- if(g=='CUP') 'S2_CUP_Adeno_2004_2022_labels.txt' else paste0('S2_',g,'_Adeno_Distant_2004_2022_CSS.txt')
 x <- fread(file.path(src,f),colClasses='character',strip.white=TRUE)
 setnames(x,c('Patient ID','Record number recode','SEER cause-specific death classification'),c('id','record','raw'))
 stopifnot(!anyDuplicated(x[,.(id,record)]))
 rec <- c('Alive or dead of other cause'=0L,'Dead (attributable to this cancer dx)'=1L,'Dead (missing/unknown COD)'=8L,'N/A not seq 0-59'=9L,'0'=0L,'1'=1L,'8'=8L,'9'=9L)
 x[,css:=unname(rec[raw])]
 e[group==g,css:=x[.SD,on=.(id,record),css],.SDcols=c('id','record')]
}
stopifnot(all(e[keep,css]%in%c(0,1,8)),all(e[keep&event==0,css]==0))
fwrite(e[keep,.(n=.N,deaths=sum(event),attributed=sum(css==1),other=sum(event==1&css==0),unknown=sum(css==8),unknown_1yr=sum(css==8&t<=365.25),deaths_1yr=sum(event==1&t<=365.25)),by=.(group,era)][order(group,era)],file.path(out,'css_coverage.csv'))
ms <- list()
for(er in eras)for(g in groups){
 lab <- paste(er,g,sep='__')
 m <- build_model(which(keep&e$group=='CUP'),which(keep&e$group==g&as.character(e$era)==er),lab,'Glandular',er,g)
 m$css <- e$css[m$source]
 ms[[lab]] <- m
}
calc <- function(m,counts=NULL){
 f<-fit_model(m,counts,diagnostic=TRUE)
 stopifnot(!nzchar(f$error))
 w<-f$odds
 if(!is.null(counts))w<-w*counts[m$source]
 a<-rowsum(cbind(w,w*(m$css==1),w*(m$css%in%c(1,8))),m$timebin,reorder=FALSE)
 risk<-rev(cumsum(rev(a[,1])))
 j<-m$times<=365.25 & risk>0
 c(os=f$values[['OS365']],css_censor=prod(1-a[j,2]/risk[j]),css_event=prod(1-a[j,3]/risk[j]))
}
point <- t(vapply(ms,calc,numeric(3)))
for(lab in names(ms)){
 m<-ms[[lab]]; w<-fit_model(m,diagnostic=TRUE)$odds
 for(k in c('css_censor','css_event')){
  ev<-if(k=='css_censor')m$css==1 else m$css%in%c(1,8)
  fit<-survival::survfit(survival::Surv(e$t[m$source],ev)~1,weights=w)
  stopifnot(abs(summary(fit,times=365.25)$surv-point[lab,k])<1e-11)
 }
}
repfit<-function(b){
 set.seed(20261025+b)
 counts<-integer(N)
 for(idx in strata) counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx))
 t(vapply(ms,calc,numeric(3),counts=counts))
}
cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS', '2')))
clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL})
clusterExport(cl,c('N','strata','ms','repfit','calc','fit_model','km_cached'))
B<-2000L; boot<-list()
for(i in seq(1L,B,by=100L)){
 boot<-c(boot,parLapplyLB(cl,i:min(i+99L,B),repfit))
 cat(length(boot),'\n');flush.console()
}
stopCluster(cl)
A<-array(unlist(boot),c(10,3,B),dimnames=list(rownames(point),colnames(point),NULL))
os_draws<-readRDS('work/glandular_models/summary_arrays.rds')$A
stopifnot(all(is.finite(A)),max(abs(t(A[,'os',])-os_draws[,paste0('Era__Glandular__',rownames(point)),'OS365']))<1e-11,all(A[,'css_event',]<=A[,'css_censor',]+1e-12))
summ<-function(p,x)c(estimate=100*p,lo=100*unname(quantile(x,.025)),hi=100*unname(quantile(x,.975)))
res<-list(); absres<-list()
for(k in colnames(point)){
 for(lab in rownames(point))absres[[length(absres)+1]]<-as.data.table(c(list(model=k,group_era=lab),as.list(summ(point[lab,k],A[lab,k,]))))
 for(g in groups[-1]){
 labs<-paste(rep(rev(eras),each=2),c('CUP',g),sep='__'); sgn<-c(1,-1,-1,1)
 res[[length(res)+1]]<-as.data.table(c(list(model=k,control=g),as.list(summ(sum(point[labs,k]*sgn),colSums(A[labs,k,]*sgn)))))
 }
}
fwrite(rbindlist(res),file.path(out,'css_gap_changes.csv'))
fwrite(rbindlist(absres),file.path(out,'css_absolute.csv'))
saveRDS(list(point=point, draws=A), 'work/cause_specific_arrays.rds')
print(rbindlist(res))

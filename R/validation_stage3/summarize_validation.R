# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'));out<-file.path(root,'outputs/validation_stage3')
pt<-fread(file.path(out,'01_fixed_target_point_estimates.csv'))
boot<-readRDS(file.path(out,'bootstrap_replicates.rds'))
models<-readRDS(file.path(root,'work/stage3_model_cache.rds'))
measures<-unique(pt$measure)
meta<-unique(pt[,.(label,definition,era,group)])
stopifnot(identical(meta$label,names(models)),length(boot)>=1000)
A<-array(NA_real_,dim=c(length(boot),nrow(meta),length(measures)),dimnames=list(NULL,meta$label,measures))
for(b in seq_along(boot))A[b,,]<-boot[[b]]$v
summary_rows<-list();gap_rows<-list();change_rows<-list();coding_rows<-list();shape_rows<-list()
ci<-function(x){v<-x[is.finite(x)];c(n=length(v),lo=quantile(v,.025,names=FALSE),hi=quantile(v,.975,names=FALSE),lo_family=quantile(v,.00625,names=FALSE),hi_family=quantile(v,.99375,names=FALSE),sd=sd(v))}
for(i in seq_len(nrow(meta)))for(k in seq_along(measures)){
 p<-pt[label==meta$label[i]&measure==measures[k],estimate]
 summary_rows[[length(summary_rows)+1]]<-cbind(meta[i],data.table(measure=measures[k],estimate=p),as.data.table(as.list(ci(A[,i,k]))))
}
est<-rbindlist(summary_rows)
for(def in c('Broad','Strict8140'))for(er in c('2004-2007','2015-2019','2020-2022'))for(g in c('Lung','Pancreas','Colorectal','Stomach')){
 ic<-which(meta$definition==def&meta$era==er&meta$group=='CUP');ig<-which(meta$definition==def&meta$era==er&meta$group==g)
 for(k in seq_along(measures)){
  gp<-pt[label==meta$label[ic]&measure==measures[k],estimate]-pt[label==meta$label[ig]&measure==measures[k],estimate]
  gb<-A[,ic,k]-A[,ig,k]
  gap_rows[[length(gap_rows)+1]]<-cbind(data.table(definition=def,era=er,control=g,measure=measures[k],estimate=gp),as.data.table(as.list(ci(gb))))
  if(er!='2004-2007'){
   ec<-which(meta$definition==def&meta$era=='2004-2007'&meta$group=='CUP');eg<-which(meta$definition==def&meta$era=='2004-2007'&meta$group==g)
   ep<-pt[label==meta$label[ec]&measure==measures[k],estimate]-pt[label==meta$label[eg]&measure==measures[k],estimate]
   bc<-gb-(A[,ec,k]-A[,eg,k]);idx<-seq_along(bc)
   half<-length(bc)/2
   change_rows[[length(change_rows)+1]]<-cbind(data.table(definition=def,era=er,control=g,measure=measures[k],estimate=gp-ep,split_CI_endpoint_max_difference=max(abs(quantile(bc[idx<=half],c(.025,.975),na.rm=TRUE)-quantile(bc[idx>half],c(.025,.975),na.rm=TRUE)))),as.data.table(as.list(ci(bc))))
  }
 }
 # Change in the absolute gap from day90 to day365 is endpoint dependence, NOT a PH test.
 shape_point<-(pt[label==meta$label[ic]&measure=='OS365',estimate]-pt[label==meta$label[ig]&measure=='OS365',estimate])-(pt[label==meta$label[ic]&measure=='OS90',estimate]-pt[label==meta$label[ig]&measure=='OS90',estimate])
 sb<-(A[,ic,'OS365']-A[,ig,'OS365'])-(A[,ic,'OS90']-A[,ig,'OS90'])
 shape_rows[[length(shape_rows)+1]]<-cbind(data.table(definition=def,era=er,control=g,contrast='OS365_gap_minus_OS90_gap',estimate=shape_point),as.data.table(as.list(ci(sb))))
}
for(er in c('2015-2019','2020-2022'))for(g in c('Lung','Pancreas','Colorectal','Stomach')){
 bc<-function(def){ic<-which(meta$definition==def&meta$era==er&meta$group=='CUP');ig<-which(meta$definition==def&meta$era==er&meta$group==g);ec<-which(meta$definition==def&meta$era=='2004-2007'&meta$group=='CUP');eg<-which(meta$definition==def&meta$era=='2004-2007'&meta$group==g);A[,ic,'OS365']-A[,ig,'OS365']-A[,ec,'OS365']+A[,eg,'OS365']}
 changes<-rbindlist(change_rows)
 p<-changes[definition=='Strict8140'&era==er&control==g&measure=='OS365',estimate]-changes[definition=='Broad'&era==er&control==g&measure=='OS365',estimate]
 coding_rows[[length(coding_rows)+1]]<-cbind(data.table(era=er,control=g,contrast='Strict_minus_broad_gap_change',estimate=p),as.data.table(as.list(ci(bc('Strict8140')-bc('Broad')))))
}
fwrite(est,file.path(out,'05_group_era_bootstrap_CI.csv'),bom=TRUE)
fwrite(rbindlist(gap_rows),file.path(out,'06_gap_bootstrap_CI.csv'),bom=TRUE)
fwrite(rbindlist(change_rows),file.path(out,'07_gap_change_bootstrap_CI.csv'),bom=TRUE)
fwrite(rbindlist(coding_rows),file.path(out,'08_coding_change_bootstrap_CI.csv'),bom=TRUE)
fwrite(rbindlist(shape_rows),file.path(out,'09_time_window_bootstrap_CI.csv'),bom=TRUE)
ch<-rbindlist(change_rows)
print(ch[definition=='Broad'&era=='2020-2022'&measure=='OS365',.(control,estimate_pp=100*estimate,lo_pp=100*lo,hi_pp=100*hi,lo_family_pp=100*lo_family,hi_family_pp=100*hi_family,n)])
writeLines('Confidence intervals include reestimated weights and shared CUP reference resampling; post-screening validation,not independent confirmation.',file.path(out,'CI_COMPLETE.txt'))

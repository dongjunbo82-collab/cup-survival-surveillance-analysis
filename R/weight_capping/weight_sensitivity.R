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
suppressPackageStartupMessages({library(data.table);library(splines);library(parallel);library(jsonlite)})
setDTthreads(1);root<-normalizePath(getwd(),winslash='/');out<-file.path(root,'outputs/weight_capping');dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_json(list(date='2026-11-02',status='Post hoc sensitivity; primary analysis unchanged',rule='Winsorize source odds at within-model 1st and 99th empirical percentiles, type 7; recompute with bootstrap multiplicities before multiplying odds by multiplicity',B=2000,seed_base=20261025,references=c('Pooled','Early')),file.path(out,'analysis_specification.json'),pretty=TRUE)
e<-readRDS('work/stage2_eligible.rds');source('R/method_reproducibility/analysis_functions.R');cfg<-fromJSON('config/histology_definitions.json');keep<-e$hist%in%cfg$reviewed_keep_codes
fit_original<-fit_model
fit_model<-function(m,counts=NULL,diagnostic=FALSE)fit_original(m,counts,diagnostic=TRUE)
weighted_q<-function(x,n,p=c(.01,.99)){j<-n>0;x<-x[j];n<-n[j];o<-order(x);x<-x[o];cs<-cumsum(n[o]);h<-(sum(n)-1)*p+1;lo<-floor(h);hi<-ceiling(h);xl<-x[findInterval(lo-1,cs)+1L];xh<-x[findInterval(hi-1,cs)+1L];xl+(h-lo)*(xh-xl)}
stopifnot(max(abs(weighted_q(c(1,4,7),c(3,2,8))-quantile(rep(c(1,4,7),c(3,2,8)),c(.01,.99),names=FALSE)))<1e-12)
models<-list();diag<-list();points<-list()
for(ref in c('Pooled','Early'))for(er in c('2004-2007','2020-2022'))for(g in c('CUP','Lung','Pancreas','Colorectal','Stomach')){
 target<-which(keep&e$group=='CUP'&(ref=='Pooled'|as.character(e$era)=='2004-2007'));lab<-paste(ref,er,g,sep='__');m<-build_model(target,which(keep&e$group==g&as.character(e$era)==er),lab,ref,er,g);models[[lab]]<-m
 f<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(f$error));q<-weighted_q(f$odds,rep(1,length(f$odds)));cw<-pmin(pmax(f$odds,q[1]),q[2]);v<-km_cached(m,cw);points[[lab]]<-rbind(data.table(label=lab,type='Original',measure=measures,estimate=f$values),data.table(label=lab,type='Capped',measure=measures,estimate=v))
 X<-m$X[,-1,drop=FALSE];nt<-m$nt;ns<-m$ns;mt<-colSums(X*nt)/sum(nt);ms<-colSums(X*ns)/sum(ns);vt<-(colSums(X^2*nt)/sum(nt)-mt^2)*sum(nt)/(sum(nt)-1);vs<-(colSums(X^2*ns)/sum(ns)-ms^2)*sum(ns)/(sum(ns)-1);wb<-numeric(m$nb);rs<-rowsum(cw,m$sourcebin,reorder=TRUE);wb[as.integer(rownames(rs))]<-rs[,1];smd<-(mt-colSums(X*wb)/sum(wb))/sqrt(pmax(0,(vt+vs)/2))
 diag[[lab]]<-data.table(reference=ref,era=er,group=g,n_source=length(cw),original_ESS=f$diagnostic['ESS'],capped_ESS=sum(cw)^2/sum(cw^2),original_max_SMD=f$diagnostic['max_SMD'],capped_max_SMD=max(abs(smd),na.rm=TRUE),lower_cut=q[1],upper_cut=q[2],original_max=max(f$odds),capped_max=max(cw))
}
pt<-rbindlist(points);fwrite(pt,file.path(out,'point_survival.csv'));fwrite(rbindlist(diag),file.path(out,'weight_diagnostics.csv'))
repfit<-function(b){set.seed(20261025+b);counts<-integer(N);for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx));v<-matrix(NA_real_,length(models),6);for(i in seq_along(models)){m<-models[[i]];f<-fit_model(m,counts);if(nzchar(f$error))next;sc<-counts[m$source];q<-weighted_q(f$odds,sc);v[i,]<-km_cached(m,pmin(pmax(f$odds,q[1]),q[2])*sc)};v}
cp<-'work/cancer_control_revision/weight_capping_bootstrap.rds';boot<-if(file.exists(cp))readRDS(cp)else list();cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS','2')),outfile='work/cancer_control_revision/weight_workers.log');clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL});clusterExport(cl,c('N','strata','models','repfit','fit_model','fit_original','km_cached','weighted_q'),envir=environment())
for(start in if(length(boot)<2000)seq(length(boot)+1L,2000,by=100)else integer()){boot<-c(boot,parLapplyLB(cl,start:min(2000,start+99),repfit));saveRDS(boot,cp);cat('Completed',length(boot),'/2000\n');flush.console()};stopCluster(cl)
res<-list();sgn<-c(1,-1,-1,1)
for(ref in c('Pooled','Early'))for(g in c('Lung','Pancreas','Colorectal','Stomach')){labs<-paste(ref,c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__');ii<-match(labs,names(models));b<-vapply(boot,function(z)sum(z[ii,2]*sgn),numeric(1));p<-sum(vapply(labs,function(l)pt[label==l&type=='Capped'&measure=='OS365',estimate],numeric(1))*sgn);dd<-rbindlist(diag)[paste(reference,era,group,sep='__')%in%labs];res[[length(res)+1]]<-data.table(reference=ref,control=g,estimate_pp=p*100,lo_pp=quantile(b,.025,na.rm=TRUE)*100,hi_pp=quantile(b,.975,na.rm=TRUE)*100,lo_family_pp=quantile(b,.00625,na.rm=TRUE)*100,hi_family_pp=quantile(b,.99375,na.rm=TRUE)*100,valid_draws=sum(is.finite(b)),max_SMD=max(dd$capped_max_SMD),min_ESS=min(dd$capped_ESS))}
res<-rbindlist(res);fwrite(res,file.path(out,'capped_gap_changes.csv'));print(res)

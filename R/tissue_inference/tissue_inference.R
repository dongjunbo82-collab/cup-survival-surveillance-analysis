# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(splines);library(parallel);library(jsonlite)})
root<-normalizePath(getwd(),winslash='/');setDTthreads(1);e<-readRDS('work/stage2_eligible.rds');source('R/method_reproducibility/analysis_functions.R',encoding='UTF-8')
out<-'outputs/tissue_inference';dir.create(out,showWarnings=FALSE);cfg<-fromJSON('config/histology_definitions.json');keep<-e$hist%in%cfg$reviewed_keep_codes&e$dc3=='Histology';target<-which(keep&e$group=='CUP');groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-c('2004-2007','2020-2022');models<-list();points<-list();diags<-list();balances<-list();knots<-list()
for(er in eras)for(g in groups){
 src<-which(keep&e$group==g&as.character(e$era)==er);lab<-paste(er,g,sep='__');m<-build_model(target,src,lab,'TissueOnly',er,g);z<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(z$error));models[[lab]]<-m
 points[[lab]]<-data.table(label=lab,era=er,group=g,measure=measures,estimate=z$values)
 w<-z$odds;diags[[lab]]<-cbind(data.table(label=lab,era=er,group=g),as.data.table(as.list(c(z$diagnostic,weight_p99=unname(quantile(w,.99)),top1pct_share=sum(sort(w,decreasing=TRUE)[seq_len(ceiling(length(w)*.01))])/sum(w)))))
 balances[[lab]]<-cbind(data.table(label=lab),z$smd);aa<-ns(c(e$age[target],e$age[src]),df=4);knots[[lab]]<-data.table(label=lab,internal_knots=paste(attr(aa,'knots'),collapse=';'),boundary_knots=paste(attr(aa,'Boundary.knots'),collapse=';'))
}
pt<-rbindlist(points);dg<-rbindlist(diags);dg[,support_pass:=n_source>=100&ESS>=200&max_SMD<.1];fwrite(pt,file.path(out,'point_estimates.csv'));fwrite(dg,file.path(out,'diagnostics.csv'));fwrite(rbindlist(balances),file.path(out,'balance.csv'));fwrite(rbindlist(knots),file.path(out,'spline_settings.csv'));saveRDS(models,'work/tissue_models.rds')
# Exact verification against the existing histology-confirmed point analysis.
old<-fread('outputs/histology_revision/25_targeted_sensitivity_estimates.csv')[sensitivity=='Tissue_only'&era%in%eras];chk<-merge(pt,old,by=c('era','group','measure'));stopifnot(max(abs(chk$estimate.x-chk$estimate.y))<1e-12)
repfit<-function(b){set.seed(20261025+b);counts<-integer(N);for(idx in strata)counts[idx]<-tabulate(sample.int(length(idx),length(idx),replace=TRUE),nbins=length(idx));v<-matrix(NA_real_,length(models),6);err<-character(length(models));for(i in seq_along(models)){z<-tryCatch(fit_model(models[[i]],counts),error=function(ex)list(values=rep(NA_real_,6),error=conditionMessage(ex)));v[i,]<-z$values;err[i]<-z$error};list(b=b,v=v,error=err)}
B<-2000L;cp<-'work/tissue_bootstrap.rds';boot<-if(file.exists(cp))readRDS(cp)else list();cl<-makeCluster(as.integer(Sys.getenv('CUP_WORKERS','2')),outfile='work/tissue_workers.log');clusterEvalQ(cl,{library(data.table);setDTthreads(1);NULL});clusterExport(cl,c('N','strata','models','repfit','fit_model','km_cached'),envir=environment())
for(start in if(length(boot)<B)seq(length(boot)+1L,B,by=100)else integer()){boot<-c(boot,parLapplyLB(cl,start:min(B,start+99L),repfit));saveRDS(boot,cp);cat('Completed',length(boot),'/',B,format(Sys.time()),'\n');flush.console()};stopCluster(cl)
fail<-rbindlist(lapply(boot,function(z)data.table(replicate=z$b,label=names(models),error=z$error)));fwrite(fail[nzchar(error)],file.path(out,'bootstrap_failures.csv'));stopifnot(!any(nzchar(fail$error)));A<-array(NA_real_,c(B,length(models),6),dimnames=list(NULL,names(models),measures));for(b in 1:B)A[b,,]<-boot[[b]]$v
oldA<-readRDS('work/histology_review/summary_arrays.rds')$A;oldch<-fread('outputs/histology_revision/17_gap_changes_CI.csv')[definition=='Reviewed'&era=='2020-2022'&measure=='OS365'];res<-list()
for(g in groups[-1]){labs<-paste(c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__');sgn<-c(1,-1,-1,1);p<-sum(vapply(labs,function(l)pt[label==l&measure=='OS365',estimate],numeric(1))*sgn);v<-as.numeric(A[,labs,'OS365']%*%sgn);ol<-paste('Era','Reviewed',c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__');ov<-as.numeric(oldA[,ol,'OS365']%*%sgn)
 res[[g]]<-data.table(control=g,estimate_pp=p*100,lo_pp=100*quantile(v,.025),hi_pp=100*quantile(v,.975),lo_adjusted_pp=100*quantile(v,.00625),hi_adjusted_pp=100*quantile(v,.99375),paired_shift_pp=100*(p-oldch[control==g,estimate]),paired_shift_lo_pp=100*quantile(v-ov,.025),paired_shift_hi_pp=100*quantile(v-ov,.975),support_pass=all(dg[label%in%labs,support_pass]),max_SMD=max(dg[label%in%labs,max_SMD]),min_ESS=min(dg[label%in%labs,ESS]),max_weight=max(dg[label%in%labs,weight_max]),max_top1_share=max(dg[label%in%labs,top1pct_share]))}
fwrite(rbindlist(res),file.path(out,'gap_changes_CI.csv'));capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'));write_json(list(B=B,seed_rule='20261025+b',target_n=length(target),all_finite=all(is.finite(A)),new_primary_family=FALSE),file.path(out,'run_QA.json'),pretty=TRUE);print(rbindlist(res));cat('DONE\n')

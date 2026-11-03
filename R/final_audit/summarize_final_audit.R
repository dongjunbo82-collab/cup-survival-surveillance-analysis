# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(data.table));setDTthreads(2)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'));out<-file.path(root,'outputs/final_audit')
wt<-function(n,x)fwrite(x,file.path(out,paste0(n,'.csv')),bom=TRUE,na='NA')
pt<-fread(file.path(out,'01_point_estimates.csv'));meta<-unique(pt[,.(label,scope,definition,period,group)])
new<-readRDS(file.path(out,'bootstrap_replicates.rds'));labs<-readRDS(file.path(out,'bootstrap_labels.rds'))
old<-readRDS(file.path(root,'outputs/validation_stage3/bootstrap_replicates.rds'));oldlabs<-names(readRDS(file.path(root,'work/stage3_model_cache.rds')))
B<-length(new);stopifnot(B==2000L,length(old)==B,identical(vapply(new,function(x)x$b,integer(1)),1:B))
measures<-unique(pt$measure);A<-array(NA_real_,c(B,nrow(meta),length(measures)),dimnames=list(NULL,meta$label,measures))
ni<-match(labs,meta$label);oi<-match(paste0('Era__',oldlabs),meta$label);stopifnot(!anyNA(ni),!anyNA(oi))
for(b in seq_len(B)){A[b,ni,]<-new[[b]]$v;A[b,oi,]<-old[[b]]$v}
stopifnot(!anyNA(A))
ci<-function(x){v<-x[is.finite(x)];hh<-seq_along(x)<=length(x)/2;c(n=length(v),lo=unname(quantile(v,.025)),hi=unname(quantile(v,.975)),lo_family=unname(quantile(v,.00625)),hi_family=unname(quantile(v,.99375)),sd=sd(v),split95_endpoint_max_difference=max(abs(quantile(x[hh],c(.025,.975),na.rm=TRUE)-quantile(x[!hh],c(.025,.975),na.rm=TRUE))),split_family_endpoint_max_difference=max(abs(quantile(x[hh],c(.00625,.99375),na.rm=TRUE)-quantile(x[!hh],c(.00625,.99375),na.rm=TRUE))))}
rows<-list()
for(i in which(meta$scope=='Era'))for(k in seq_along(measures))rows[[length(rows)+1]]<-cbind(meta[i],data.table(measure=measures[k],estimate=pt[label==meta$label[i]&measure==measures[k],estimate]),as.data.table(as.list(ci(A[,i,k]))))
era_est<-rbindlist(rows);wt('22_four_era_estimates_CI',era_est)
getp<-function(lab)pt[match(paste(lab,measures),paste(label,measure)),estimate]
eras<-c('2004-2007','2008-2014','2015-2019','2020-2022');groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach')
gaps<-list();changes<-list()
for(def in c('Broad','Strict8140'))for(er in eras)for(g in groups[-1]){
 lc<-paste('Era',def,er,'CUP',sep='__');lg<-paste('Era',def,er,g,sep='__');gp<-getp(lc)-getp(lg);gb<-A[,lc,]-A[,lg,]
 for(k in seq_along(measures)){
  gaps[[length(gaps)+1]]<-cbind(data.table(definition=def,era=er,control=g,measure=measures[k],estimate=gp[k]),as.data.table(as.list(ci(gb[,k]))))
  if(er!=eras[1]){
   ec<-paste('Era',def,eras[1],'CUP',sep='__');eg<-paste('Era',def,eras[1],g,sep='__')
   cp<-gp-(getp(ec)-getp(eg));cb<-gb-(A[,ec,]-A[,eg,])
   changes[[length(changes)+1]]<-cbind(data.table(definition=def,era=er,control=g,measure=measures[k],estimate=cp[k]),as.data.table(as.list(ci(cb[,k]))))
  }
 }
}
gap<-rbindlist(gaps);change<-rbindlist(changes);wt('23_four_era_gaps_CI',gap);wt('24_four_era_gap_changes_CI',change)
shape<-list()
for(def in c('Broad','Strict8140'))for(er in eras)for(g in groups[-1]){
 lc<-paste('Era',def,er,'CUP',sep='__');lg<-paste('Era',def,er,g,sep='__')
 v<-(A[,lc,'OS365']-A[,lg,'OS365'])-(A[,lc,'OS90']-A[,lg,'OS90'])
 cp<-getp(lc)-getp(lg)
 shape[[length(shape)+1]]<-cbind(data.table(definition=def,era=er,control=g,contrast='OS365_gap_minus_OS90_gap',estimate=cp[2]-cp[1]),as.data.table(as.list(ci(v))))
}
wt('33_time_window_contrast_CI',rbindlist(shape))
# The same yearly CUP shares are used for every cancer within a definition/era.
# Conditional survival is formed from the mixture survival, NOT an average of ratios.
e<-readRDS(file.path(root,'work/stage2_eligible.rds'));dg<-fread(file.path(out,'02_model_diagnostics.csv'))
mixpt<-list();mixboot<-array(NA_real_,c(B,40L,length(measures)))
idx<-0L
for(def in c('Broad','Strict8140'))for(er in eras)for(g in groups){
 idx<-idx+1L;yrs<-sort(unique(e[as.character(era)==er,year]));labs_y<-paste('Annual',def,yrs,g,sep='__')
 q<-as.numeric(table(factor(e[group=='CUP'&as.character(era)==er&(def=='Broad'|hist==8140),year],levels=yrs)));q<-q/sum(q)
 pp<-sapply(labs_y,getp);mp<-as.numeric(pp%*%q);mp[6]<-mp[2]/mp[1]
 for(b in seq_len(B)){
  qb<-(if(def=='Broad')new[[b]]$q else new[[b]]$qs)[match(yrs,2004:2022)];qb<-qb/sum(qb)
  z<-t(A[b,labs_y,,drop=FALSE][1,,]);v<-as.numeric(z%*%qb);v[6]<-v[2]/v[1];mixboot[b,idx,]<-v
 }
 pass<-all(dg[match(labs_y,label),support_pass]);bad<-paste(yrs[!dg[match(labs_y,label),support_pass]],collapse=';')
 mixpt[[idx]]<-data.table(index=idx,definition=def,era=er,group=g,support_pass=pass,unsupported_years=bad,measure=measures,estimate=mp)
}
mix<-rbindlist(mixpt);mixci<-list()
for(i in 1:40)for(k in seq_along(measures))mixci[[length(mixci)+1]]<-cbind(mix[index==i&measure==measures[k]],as.data.table(as.list(ci(mixboot[,i,k]))))
wt('25_year_balanced_estimates_CI',rbindlist(mixci));mg<-list();mc<-list()
for(def in c('Broad','Strict8140'))for(er in eras)for(g in groups[-1]){
 ic<-unique(mix[definition==def&era==er&group=='CUP',index]);ig<-unique(mix[definition==def&era==er&group==g,index])
 p<-mix[index==ic,estimate]-mix[index==ig,estimate];v<-mixboot[,ic,]-mixboot[,ig,]
 ok<-all(mix[index%in%c(ic,ig),support_pass]);bad<-paste(unique(mix[index%in%c(ic,ig)&!support_pass,unsupported_years]),collapse=';')
 for(k in seq_along(measures)){
  mg[[length(mg)+1]]<-cbind(data.table(definition=def,era=er,control=g,measure=measures[k],estimate=p[k],support_pass=ok,unsupported_years=bad),as.data.table(as.list(ci(v[,k]))))
  if(er!=eras[1]){
   ec<-unique(mix[definition==def&era==eras[1]&group=='CUP',index]);eg<-unique(mix[definition==def&era==eras[1]&group==g,index]);pv<-p-(mix[index==ec,estimate]-mix[index==eg,estimate]);bv<-v-(mixboot[,ec,]-mixboot[,eg,])
   valid<-ok&&all(mix[index%in%c(ec,eg),support_pass])
   mc[[length(mc)+1]]<-cbind(data.table(definition=def,era=er,control=g,measure=measures[k],estimate=pv[k],support_pass=valid),as.data.table(as.list(ci(bv[,k]))))
  }
 }
}
mg<-rbindlist(mg);mc<-rbindlist(mc);wt('26_year_balanced_gaps_CI',mg);wt('27_year_balanced_gap_changes_CI',mc)
cmp<-merge(change[measure=='OS365'],mc[measure=='OS365'],by=c('definition','era','control','measure'),suffixes=c('_era','_yearbalanced'));cmp[,year_mix_shift_pp:=100*(estimate_yearbalanced-estimate_era)]
wt('28_year_mix_effect_on_primary',cmp)
fails<-fread(file.path(out,'13_bootstrap_failures.csv'));success<-mean(vapply(new,function(z)all(!nzchar(z$error)),logical(1)))
stopifnot(success>=.95,nrow(fails)==0L)
saveRDS(list(meta=meta,A=A,year_mix=mix,year_mix_boot=mixboot),file.path(root,'work/final_audit_summary_arrays.rds'))
capture.output(print(cmp[definition=='Broad'&era=='2020-2022',.(control,estimate_era,lo_era,hi_era,estimate_yearbalanced,lo_yearbalanced,hi_yearbalanced,year_mix_shift_pp,support_pass)]),file=file.path(out,'PRIMARY_RESULTS_console.txt'))
writeLines(paste('Bootstrap complete;',B,'replicates;',length(labs),'additional models per replicate; all-model success',success),file.path(out,'CI_COMPLETE.txt'))
print(cmp[definition=='Broad'&era=='2020-2022',.(control,estimate_era,estimate_yearbalanced,year_mix_shift_pp,support_pass)])

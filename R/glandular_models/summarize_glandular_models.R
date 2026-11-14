# Generated artifacts are written to local output directories.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/glandular_models', 'work/landmark',
            'outputs/standardization', 'outputs/cohort_checks',
            'outputs/definition_models', 'outputs/era_models', 'outputs/glandular_models',
            'outputs/pandemic_sensitivity', 'outputs/absolute_outcomes',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
library(data.table);setDTthreads(1)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'));out<-file.path(root,'outputs/glandular_models')
wt<-function(n,x)fwrite(x,file.path(out,paste0(n,'.csv')),bom=TRUE,na='NA')
pt<-fread(file.path(out,'04_point_estimates.csv'));meta<-unique(pt[,.(label,scope,definition,period,group)])
boot<-readRDS(file.path(out,'bootstrap_replicates.rds'));labs<-readRDS(file.path(out,'bootstrap_labels.rds'));B<-length(boot);stopifnot(B==2000,identical(labs,meta$label))
measures<-unique(pt$measure);A<-array(NA_real_,c(B,length(labs),6),dimnames=list(NULL,labs,measures));for(b in 1:B)A[b,,]<-boot[[b]]$v
stopifnot(!anyNA(A));ci<-function(x)c(n=sum(is.finite(x)),lo=unname(quantile(x,.025)),hi=unname(quantile(x,.975)),lo_family=unname(quantile(x,.00625)),hi_family=unname(quantile(x,.99375)),split95=max(abs(quantile(x[1:1000],c(.025,.975))-quantile(x[1001:2000],c(.025,.975)))))
getp<-function(lab)pt[match(paste(lab,measures),paste(label,measure)),estimate]
rows<-list();gaps<-list();changes<-list();shape<-list();eras<-c('2004-2007','2008-2014','2015-2019','2020-2022');groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach')
for(i in which(meta$scope=='Era'))for(k in 1:6)rows[[length(rows)+1]]<-cbind(meta[i],data.table(measure=measures[k],estimate=getp(labs[i])[k]),as.data.table(as.list(ci(A[,i,k]))))
wt('15_era_estimates_CI',rbindlist(rows))
for(def in c('Glandular','MinimalNoNE'))for(er in eras)for(g in groups[-1]){
 lc<-paste('Era',def,er,'CUP',sep='__');lg<-paste('Era',def,er,g,sep='__');p<-getp(lc)-getp(lg);v<-A[,lc,]-A[,lg,]
 for(k in 1:6){
  gaps[[length(gaps)+1]]<-cbind(data.table(definition=def,era=er,control=g,measure=measures[k],estimate=p[k]),as.data.table(as.list(ci(v[,k]))))
  if(er!=eras[1]){
   ec<-paste('Era',def,eras[1],'CUP',sep='__');eg<-paste('Era',def,eras[1],g,sep='__');pv<-p-(getp(ec)-getp(eg));bv<-v-(A[,ec,]-A[,eg,])
   changes[[length(changes)+1]]<-cbind(data.table(definition=def,era=er,control=g,measure=measures[k],estimate=pv[k]),as.data.table(as.list(ci(bv[,k]))))
  }
 }
 shape[[length(shape)+1]]<-cbind(data.table(definition=def,era=er,control=g,estimate=p[2]-p[1]),as.data.table(as.list(ci(v[,'OS365']-v[,'OS90']))))
}
wt('16_gaps_CI',rbindlist(gaps));ch<-rbindlist(changes);wt('17_gap_changes_CI',ch);wt('18_time_window_contrasts_CI',rbindlist(shape))
e<-readRDS(file.path(root,'work/eligible_cohort.rds'));cfg<-jsonlite::fromJSON(file.path(root,'config/histology_definitions.json'));dg<-fread(file.path(out,'05_weight_diagnostics.csv'))
mix<-list();mb<-array(NA_real_,c(B,20,6));idx<-0L
for(er in eras)for(g in groups){
 idx<-idx+1L;yrs<-sort(unique(e[as.character(era)==er,year]));yl<-paste('Annual','Glandular',yrs,g,sep='__');q<-as.numeric(table(factor(e[group=='CUP'&hist%in%cfg$glandular_codes&as.character(era)==er,year],levels=yrs)));q<-q/sum(q)
 p<-as.numeric(sapply(yl,getp)%*%q);p[6]<-p[2]/p[1]
 for(b in 1:B){qb<-boot[[b]]$q[match(yrs,2004:2022)];qb<-qb/sum(qb);v<-as.numeric(t(A[b,yl,])%*%qb);v[6]<-v[2]/v[1];mb[b,idx,]<-v}
 mix[[idx]]<-data.table(index=idx,era=er,group=g,measure=measures,estimate=p,support_pass=all(dg[match(yl,label),support_pass]),unsupported_years=paste(yrs[!dg[match(yl,label),support_pass]],collapse=';'))
}
mix<-rbindlist(mix);mixci<-list();mg<-list();mc<-list()
for(i in 1:20)for(k in 1:6)mixci[[length(mixci)+1]]<-cbind(mix[index==i&measure==measures[k]],as.data.table(as.list(ci(mb[,i,k]))))
wt('19_year_balanced_estimates_CI',rbindlist(mixci))
for(er in eras)for(g in groups[-1]){
 ic<-unique(mix[era==er&group=='CUP',index]);ig<-unique(mix[era==er&group==g,index]);gp<-mix[index==ic,estimate]-mix[index==ig,estimate];gb<-mb[,ic,]-mb[,ig,]
 if(er==eras[1])next
 ec<-unique(mix[era==eras[1]&group=='CUP',index]);eg<-unique(mix[era==eras[1]&group==g,index]);p<-gp-mix[index==ec,estimate]+mix[index==eg,estimate];v<-gb-mb[,ec,]+mb[,eg,]
 for(k in 1:6)mc[[length(mc)+1]]<-cbind(data.table(era=er,control=g,measure=measures[k],estimate=p[k],support_pass=all(mix[index%in%c(ic,ig,ec,eg),support_pass])),as.data.table(as.list(ci(v[,k]))))
}
wt('20_year_balanced_gap_changes_CI',rbindlist(mc))
# Paired changes relative to the broad morphology definition.
reference_result<-readRDS(file.path(root,'outputs/definition_models/bootstrap_replicates.rds'));ol<-names(readRDS(file.path(root,'work/definition_model_cache.rds')));oldch<-fread(file.path(root,'outputs/definition_models/07_gap_change_bootstrap_CI.csv'))
comparison<-list()
for(def in c('Glandular','MinimalNoNE'))for(g in groups[-1]){
 newlabs<-paste('Era',def,c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__');oldlabs<-paste('Broad',c('2020-2022','2020-2022','2004-2007','2004-2007'),c('CUP',g,'CUP',g),sep='__')
 bv<-A[,newlabs[1],'OS365']-A[,newlabs[2],'OS365']-A[,newlabs[3],'OS365']+A[,newlabs[4],'OS365']
 ov<-vapply(reference_result,function(z)sum(z$v[match(oldlabs,ol),2]*c(1,-1,-1,1)),numeric(1));np<-ch[definition==def&era=='2020-2022'&control==g&measure=='OS365',estimate];op<-oldch[definition=='Broad'&era=='2020-2022'&control==g&measure=='OS365',estimate]
 comparison[[length(comparison)+1]]<-cbind(data.table(definition=def,control=g,broad=op,glandular=np,estimate=np-op),as.data.table(as.list(ci(bv-ov))))
}
wt('21_paired_effect_of_definition_comparison',rbindlist(comparison));wt('22_strict8140_reference',oldch[definition=='Strict8140'])
saveRDS(list(A=A,meta=meta,mix=mix,mb=mb),file.path(root,'work/glandular_models/summary_arrays.rds'))
print(ch[definition=='Glandular'&era=='2020-2022'&measure=='OS365']);writeLines('All glandular intervals complete',file.path(out,'CI_COMPLETE.txt'))

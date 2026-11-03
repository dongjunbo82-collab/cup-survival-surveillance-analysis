# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(splines)})
setDTthreads(1);out<-'outputs/early_landmark';eras<-c('2004-2007','2008-2014','2015-2019','2020-2022');pt<-fread('outputs/histology_revision/15_era_estimates_CI.csv')[definition=='Reviewed'&group%in%c('CUP','Pancreas')&measure%in%c('OS90','OS365','Conditional365given90')];pt[measure=='OS90',`:=`(measure='Death90',estimate=1-estimate,lo=1-hi,hi=1-lo)];fwrite(pt[,.(label,era=period,group,measure,estimate,lo,hi,draws=n)],file.path(out,'four_era_absolute_probabilities.csv'))
cp<-'work/cancer_control_revision/standalone_reproduction.rds';stopifnot(file.exists(cp));boot<-readRDS(cp);stopifnot(length(boot)==2000);res<-list();vlist<-list();get<-function(g,er,me){i<-match(g,c('CUP','Lung','Pancreas','Colorectal','Stomach'))+if(er==eras[4])5L else 0L;j<-match(if(me=='Death90')'OS90'else me,c('OS90','OS365','RMST90','RMST365','RMST90to365','Conditional365given90'));v<-vapply(boot,function(x)x[i,j],numeric(1));if(me=='Death90')1-v else v}
for(g in c('CUP','Pancreas'))for(me in c('Death90','OS365','Conditional365given90')){v<-get(g,eras[4],me)-get(g,eras[1],me);p<-pt[group==g&period==eras[4]&measure==me,estimate]-pt[group==g&period==eras[1]&measure==me,estimate];vlist[[paste(g,me)]]<-v;res[[length(res)+1]]<-data.table(comparison=g,measure=me,estimate_pp=p*100,lo_pp=quantile(v,.025)*100,hi_pp=quantile(v,.975)*100,draws=length(v))}
for(me in c('Death90','OS365','Conditional365given90')){v<-vlist[[paste('CUP',me)]]-vlist[[paste('Pancreas',me)]];rr<-rbindlist(res);p<-rr[comparison=='CUP'&measure==me,estimate_pp]-rr[comparison=='Pancreas'&measure==me,estimate_pp];res[[length(res)+1]]<-data.table(comparison='CUP-minus-pancreas change contrast',measure=me,estimate_pp=p,lo_pp=quantile(v,.025)*100,hi_pp=quantile(v,.975)*100,draws=length(v))}
fwrite(rbindlist(res),file.path(out,'first_last_era_changes.csv'))
# Diagnose outcome heterogeneity in the high-weight tail without replacing primary curves.
e<-readRDS('work/stage2_eligible.rds');source('R/method_reproducibility/analysis_functions.R');cfg<-fromJSON('config/histology_definitions.json');keep<-e$hist%in%cfg$reviewed_keep_codes;ans<-list()
for(ref in c('Pooled','Early'))for(er in eras[c(1,4)])for(g in c('CUP','Lung','Pancreas','Colorectal','Stomach')){tt<-which(keep&e$group=='CUP'&(ref=='Pooled'|as.character(e$era)==eras[1]));m<-build_model(tt,which(keep&e$group==g&as.character(e$era)==er),'diag',ref,er,g);f<-fit_model(m,diagnostic=TRUE);w<-f$odds;n<-length(w);hi<-order(w,decreasing=TRUE)[seq_len(ceiling(.01*n))];for(sub in c('All','Top 1%','Other 99%')){sel<-if(sub=='All')rep(TRUE,n)else if(sub=='Top 1%')seq_len(n)%in%hi else !seq_len(n)%in%hi;v<-km_cached(m,w*sel);ans[[length(ans)+1]]<-data.table(reference=ref,era=er,group=g,subset=sub,source_n=sum(sel),weight_fraction=sum(w[sel])/sum(w),OS90=v[1],OS365=v[2])}}
fwrite(rbindlist(ans),file.path(out,'weight_tail_outcomes.csv'));print(rbindlist(res))

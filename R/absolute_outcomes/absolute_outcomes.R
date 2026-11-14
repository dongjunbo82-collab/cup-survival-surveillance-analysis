# Generated artifacts are written to local output directories.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/glandular_models', 'work/landmark',
            'outputs/standardization', 'outputs/cohort_checks',
            'outputs/definition_models', 'outputs/era_models', 'outputs/glandular_models',
            'outputs/pandemic_sensitivity', 'outputs/absolute_outcomes',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
library(data.table);library(jsonlite);setDTthreads(1)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'))
out<-file.path(root,'outputs/absolute_outcomes');dir.create(out,showWarnings=FALSE)
src<-file.path(root,'outputs/glandular_models');pt<-fread(file.path(src,'04_point_estimates.csv'))
A<-readRDS(file.path(root,'work/glandular_models/summary_arrays.rds'))$A
groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-c('2004-2007','2008-2014','2015-2019','2020-2022')
lab<-function(er,g)paste('Era','Glandular',er,g,sep='__')
getp<-function(l,m)pt[label==l&measure==m,estimate]
ci<-function(v)c(lo=unname(quantile(v,.025)),hi=unname(quantile(v,.975)))
rows<-list();changes<-list()
for(g in groups){
 a<-lab(eras[1],g);b<-lab(eras[4],g)
 for(er in eras){l<-lab(er,g);rows[[length(rows)+1]]<-cbind(data.table(group=g,era=er,estimate=getp(l,'OS365')),as.data.table(as.list(ci(A[,l,'OS365']))))}
 changes[[g]]<-cbind(data.table(group=g,estimate=getp(b,'OS365')-getp(a,'OS365')),as.data.table(as.list(ci(A[,b,'OS365']-A[,a,'OS365']))))
}
abs<-rbindlist(rows);chg<-rbindlist(changes);fwrite(abs,file.path(out,'36_absolute_OS365.csv'));fwrite(chg,file.path(out,'37_within_group_OS365_changes.csv'))
latest<-list()
for(g in c('CUP','Pancreas'))for(m in c('Death90','OS365','Conditional365given90')){
 l<-lab(eras[4],g);mm<-if(m=='Death90')'OS90'else m
 v<-A[,l,mm];p<-getp(l,mm);if(m=='Death90'){v<-1-v;p<-1-p}
 latest[[length(latest)+1]]<-cbind(data.table(group=g,measure=m,estimate=p),as.data.table(as.list(ci(v))))
}
latest<-rbindlist(latest);fwrite(latest,file.path(out,'38_latest_absolute_probabilities.csv'))
e<-as.data.table(readRDS(file.path(root,'work/eligible_cohort.rds')));cfg<-fromJSON(file.path(root,'config/histology_definitions.json'));e<-e[hist%in%cfg$glandular_codes]
cup<-e[group=='CUP'];cuprows<-list()
for(er in eras){z<-cup[as.character(era)==er];n<-nrow(z)
 vals<-list(n=n,age_median=median(z$age),age_q1=unname(quantile(z$age,.25)),age_q3=unname(quantile(z$age,.75)),age75=sum(z$age>=75),female=sum(z$sex=='Female'),histology=sum(z$dc3=='Histology'),cytology=sum(z$dc3=='Cytology'),micro_unspecified=sum(z$dc3=='Micro_unspecified'),strict8140=sum(z$hist==8140),marital_unknown=sum(z$marital3=='Unknown'),race_unknown=sum(grepl('Unknown',z$race,ignore.case=TRUE)|is.na(z$race)|z$race==''),age_unknown=sum(is.na(z$age)),sex_unknown=sum(is.na(z$sex)|!z$sex%in%c('Male','Female')),income_unknown=sum(!is.finite(z$income_num)),rural_unknown=sum(z$rural3=='Unknown'|is.na(z$rural3)))
 for(v in names(vals))cuprows[[length(cuprows)+1]]<-data.table(era=er,variable=v,value=vals[[v]],denominator=n,percent=if(v%in%c('n','age_median','age_q1','age_q3'))NA_real_ else 100*vals[[v]]/n)
}
fwrite(rbindlist(cuprows),file.path(out,'39_CUP_era_composition.csv'))
fwrite(e[,.(n=.N,followup_year_min=min(as.integer(followyear),na.rm=TRUE),followup_year_max=max(as.integer(followyear),na.rm=TRUE),diagnosis_min=min(year),diagnosis_max=max(year)),by=group],file.path(out,'40_followup_and_year_QA.csv'))
stopifnot(nrow(e)==380648,nrow(cup)==27715,all(dim(A)==c(2000,135,6))||dim(A)[1]==2000)
print(latest);print(chg);print(rbindlist(cuprows)[variable%in%c('age_median','histology','strict8140','marital_unknown','race_unknown')])

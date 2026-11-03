# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
Sys.setenv(STAGE3_FUNCTIONS_ONLY='1');source('R/validation_stage3/run_validation.R',encoding='UTF-8');Sys.unsetenv('STAGE3_FUNCTIONS_ONLY')
library(jsonlite);setDTthreads(1);out<-file.path(root,'outputs/histology_revision');cfg<-fromJSON(file.path(root,'config/histology_definitions.json'));base<-e$hist%in%cfg$reviewed_keep_codes;base_t<-e$t
rr<-list();dd<-list()
for(sens in c('Tissue_only','Complete_dates','Month_scale','Exclude_appendix','Socioeconomic_expanded')){
 use<-base&switch(sens,Tissue_only=e$dc3=='Histology',Complete_dates=grepl('^Complete',e$flag),Exclude_appendix=!(e$group=='Colorectal'&startsWith(e$site,'C18.1-')),rep(TRUE,N))
 e[,t:=base_t];if(sens=='Month_scale')e[,t:=pmax(months,.5)*365.25/12]
 target<-which(use&e$group=='CUP')
 for(er in levels(e$era))for(g in c('CUP','Lung','Pancreas','Colorectal','Stomach')){
  source<-which(use&e$group==g&as.character(e$era)==er);lab<-paste(sens,er,g,sep='__');m<-build_model(target,source,lab,sens,er,g)
  if(sens=='Socioeconomic_expanded'){
   zz<-rbind(e[target],e[m$source]);xx<-model.matrix(~ns(age,df=4)+old90+factor(sex)+factor(race)+factor(marital3)+factor(dc3)+factor(income4)+factor(rural3),zz)
   keys<-do.call(paste,c(zz[,c(covnames,'income4','rural3'),with=FALSE],sep='|'));uk<-unique(keys);bins<-match(keys,uk)
   m$X<-xx[match(uk,keys),,drop=FALSE];m$nb<-length(uk);m$targetbin<-head(bins,length(target));m$sourcebin<-tail(bins,length(source));m$nt<-tabulate(m$targetbin,m$nb);m$ns<-tabulate(m$sourcebin,m$nb)
  }
  f<-fit_model(m,diagnostic=TRUE);stopifnot(!nzchar(f$error));rr[[lab]]<-data.table(sensitivity=sens,era=er,group=g,measure=measures,estimate=f$values);dd[[lab]]<-cbind(data.table(sensitivity=sens,era=er,group=g),as.data.table(as.list(f$diagnostic)))
 }
}
r<-rbindlist(rr);d<-rbindlist(dd);d[,support_pass:=ESS>=200&max_SMD<.1]
cup<-r[group=='CUP',.(sensitivity,era,measure,CUP=estimate)];gaps<-merge(r[group!='CUP'],cup,by=c('sensitivity','era','measure'));gaps[,gap:=CUP-estimate]
early<-gaps[era=='2004-2007',.(sensitivity,group,measure,early=gap)];ch<-merge(gaps[era!='2004-2007'],early,by=c('sensitivity','group','measure'));ch[,gap_change:=gap-early]
fwrite(r,file.path(out,'25_targeted_sensitivity_estimates.csv'),bom=TRUE);fwrite(d,file.path(out,'26_targeted_sensitivity_diagnostics.csv'),bom=TRUE);fwrite(ch,file.path(out,'27_targeted_sensitivity_changes.csv'),bom=TRUE)
z<-e[base];tab<-z[,.(n=.N,age_median=median(age),age_q1=quantile(age,.25),age_q3=quantile(age,.75),female=sum(sex=='Female'),older75=sum(age>=75),histology=sum(dc3=='Histology'),cytology=sum(dc3=='Cytology'),strict8140=sum(hist==8140),alive_censored_under1y=sum(event==0&days<365.25)),by=group]
fwrite(tab,file.path(out,'28_baseline_summary.csv'),bom=TRUE)
fwrite(z[,.(n=.N,zero_days=sum(days==0),zero_months=sum(months==0),alive_censor1y=sum(event==0&days<365.25),histology_pct=100*mean(dc3=='Histology'),strict8140_pct=100*mean(hist==8140)),by=.(group,era)],file.path(out,'29_cohort_era_description.csv'),bom=TRUE)
print(ch[era=='2020-2022'&measure=='OS365',.(sensitivity,group,gap_change)])

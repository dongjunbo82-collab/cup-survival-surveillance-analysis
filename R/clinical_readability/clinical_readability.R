# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
library(data.table);library(jsonlite);setDTthreads(1)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'))
out<-file.path(root,'outputs/clinical_readability');dir.create(out,showWarnings=FALSE)
src<-file.path(root,'outputs/histology_revision');pt<-fread(file.path(src,'04_point_estimates.csv'))
A<-readRDS(file.path(root,'work/histology_review/summary_arrays.rds'))$A
groups<-c('CUP','Lung','Pancreas','Colorectal','Stomach');eras<-c('2004-2007','2008-2014','2015-2019','2020-2022')
lab<-function(er,g)paste('Era','Reviewed',er,g,sep='__')
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
e<-as.data.table(readRDS(file.path(root,'work/stage2_eligible.rds')));cfg<-fromJSON(file.path(root,'config/histology_definitions.json'));e<-e[hist%in%cfg$reviewed_keep_codes]
cup<-e[group=='CUP'];cuprows<-list()
for(er in eras){z<-cup[as.character(era)==er];n<-nrow(z)
 vals<-list(n=n,age_median=median(z$age),age_q1=unname(quantile(z$age,.25)),age_q3=unname(quantile(z$age,.75)),age75=sum(z$age>=75),female=sum(z$sex=='Female'),histology=sum(z$dc3=='Histology'),cytology=sum(z$dc3=='Cytology'),micro_unspecified=sum(z$dc3=='Micro_unspecified'),strict8140=sum(z$hist==8140),marital_unknown=sum(z$marital3=='Unknown'),race_unknown=sum(grepl('Unknown',z$race,ignore.case=TRUE)|is.na(z$race)|z$race==''),age_unknown=sum(is.na(z$age)),sex_unknown=sum(is.na(z$sex)|!z$sex%in%c('Male','Female')),income_unknown=sum(!is.finite(z$income_num)),rural_unknown=sum(z$rural3=='Unknown'|is.na(z$rural3)))
 for(v in names(vals))cuprows[[length(cuprows)+1]]<-data.table(era=er,variable=v,value=vals[[v]],denominator=n,percent=if(v%in%c('n','age_median','age_q1','age_q3'))NA_real_ else 100*vals[[v]]/n)
}
fwrite(rbindlist(cuprows),file.path(out,'39_CUP_era_composition.csv'))
fwrite(e[,.(n=.N,followup_year_min=min(as.integer(followyear),na.rm=TRUE),followup_year_max=max(as.integer(followyear),na.rm=TRUE),diagnosis_min=min(year),diagnosis_max=max(year)),by=group],file.path(out,'40_followup_and_year_QA.csv'))
# New Figure 1: own survival changes and absolute early mortality/conditional probabilities.
png(file.path(out,'Figure1_absolute_survival.png'),width=2700,height=1500,res=300)
par(mfrow=c(1,2),mar=c(4.5,7.2,3,1.0),oma=c(.5,0,.5,0),family='sans',cex=.94)
cols<-c('#202020','#2166AC','#D97706','#228833','#8856A7')
plot(NA,xlim=c(0,60),ylim=c(.5,5.6),xlab='Standardized 1-year overall survival (%)',ylab='',yaxt='n');axis(2,at=5:1,labels=groups,las=1,tick=FALSE);abline(v=seq(0,60,10),col='#EEEEEE')
for(i in 1:5){g<-groups[i];a<-abs[group==g&era=='2004-2007'];b<-abs[group==g&era=='2020-2022'];y<-6-i;segments(a$estimate*100,y,b$estimate*100,y,lwd=2,col=cols[i]);for(k in 1:2){z<-if(k==1)a else b;segments(100*z$lo,y,100*z$hi,y,col=cols[i],lwd=1.5);points(z$estimate*100,y,pch=if(k==1)1 else 16,cex=1.1,col=cols[i])}}
legend('bottomright',c('2004-2007','2020-2022'),pch=c(1,16),bty='n',cex=.85);title(main='A  Survival within each cancer group',adj=0,cex.main=1)
par(mar=c(4.5,7.6,3,1))
plot(NA,xlim=c(0,85),ylim=c(.4,3.65),xlab='Standardized probability (%)',ylab='',yaxt='n');axis(2,at=3:1,labels=c('Death by 90 days','Survival at 1 year','Survival to 1 year\namong 90-day survivors'),las=1,tick=FALSE,cex.axis=.83);abline(v=seq(0,80,20),col='#EEEEEE')
for(i in 1:3)for(k in 1:2){z<-latest[measure==c('Death90','OS365','Conditional365given90')[i]&group==c('CUP','Pancreas')[k]];y<-4-i+c(.1,-.1)[k];segments(100*z$lo,y,100*z$hi,y,col=cols[c(1,3)[k]],lwd=2);points(100*z$estimate,y,col=cols[c(1,3)[k]],pch=c(16,17)[k],cex=1.1)}
legend('topright',c('CUP','Pancreas'),col=cols[c(1,3)],pch=c(16,17),bty='n',cex=.9);title(main='B  Absolute prognosis in 2020-2022',adj=0,cex.main=1)
dev.off()
stopifnot(nrow(e)==380648,nrow(cup)==27715,all(dim(A)==c(2000,135,6))||dim(A)[1]==2000)
print(latest);print(chg);print(rbindlist(cuprows)[variable%in%c('age_median','histology','strict8140','marital_unknown','race_unknown')])

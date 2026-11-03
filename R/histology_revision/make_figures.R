# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
library(data.table)
root <- Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'))
out <- file.path(root,'outputs/histology_revision')
dest <- file.path(root,'outputs/concise_manuscript')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
cur <- fread(file.path(out,'07_reviewed_curves.csv'))
risk <- fread(file.path(out,'08_risk_support.csv'))
groups <- c('CUP','Lung','Pancreas','Colorectal','Stomach')
cols <- c('#202020','#2166AC','#D97706','#228833','#8856A7')
eras <- c('2004-2007','2008-2014','2015-2019','2020-2022')
png(file.path(dest,'Figure1_survival_curves.png'),width=2400,height=2300,res=300)
layout(matrix(c(1,3,2,4,5,7,6,8),4,2,byrow=TRUE),heights=c(3.2,1.25,3.2,1.25))
par(oma=c(.7,.5,1.3,.2),family='sans',cex=.88)
for(j in seq_along(eras)){
 er<-eras[j]
 par(mar=c(3.5,4.4,2.6,1))
 plot(NA,xlim=c(0,365.25),ylim=c(0,1),xlab='',ylab='Standardized overall survival',xaxt='n',yaxs='i')
 mtext('Days since diagnosis',side=1,line=2.0,cex=.73)
 axis(1,at=c(0,90,180,270,365),labels=c(0,90,180,270,365))
 abline(h=seq(0,1,.2),col='#EEEEEE',lwd=.7)
 for(k in seq_along(groups)){z<-cur[era==er&group==groups[k]];lines(z$day,z$survival,type='s',col=cols[k],lwd=2)}
 title(main=paste(LETTERS[j],er),adj=0,cex.main=1.08)
 if(j==1)legend('topright',groups,col=cols,lty=1,lwd=2,bty='n',cex=.78)
 par(mar=c(.5,4.4,1.1,1))
 plot.new();plot.window(xlim=c(0,365.25),ylim=c(.5,6.4),xaxs='i')
 text(0,6.2,'Unweighted number at risk',adj=0,cex=.74)
 for(k in seq_along(groups)){
  text(-25,6-k,groups[k],adj=1,xpd=NA,col=cols[k],cex=.70)
  for(tau in c(0,90,365.25))text(tau,6-k,format(risk[era==er&group==groups[k]&day==tau,n_risk],big.mark=','),cex=.7,adj=c(if(tau==0)0 else if(tau>365)1 else .5,.5))
 }
}
dev.off()
gaps<-fread(file.path(out,'16_gaps_CI.csv'))[definition=='Reviewed'&era=='2020-2022'&measure%in%c('OS90','OS365')]
png(file.path(dest,'Figure2_time_window_gaps.png'),width=2400,height=1450,res=300)
par(mar=c(4.6,7,2.6,1),family='sans',cex=1)
plot(NA,xlim=c(-48,5),ylim=c(.45,4.65),xlab='CUP minus known-primary survival (percentage points)',ylab='',yaxt='n',xaxs='i')
axis(2,at=4:1,labels=groups[-1],las=1,tick=FALSE)
abline(v=0,lty=2,col='#888888');abline(v=seq(-40,-10,10),col='#EEEEEE')
for(i in 1:4)for(k in 1:2){z<-gaps[control==groups[i+1]&measure==c('OS90','OS365')[k]];y<-5-i+c(.12,-.12)[k];segments(100*z$lo,y,100*z$hi,y,col=c('#2166AC','#D97706')[k],lwd=2);points(100*z$estimate,y,pch=c(16,17)[k],cex=1.1,col=c('#2166AC','#D97706')[k])}
legend('bottomleft',c('90-day survival','1-year survival'),pch=c(16,17),col=c('#2166AC','#D97706'),bty='n',cex=.9)
title(main='2020-2022: survival differences depend on time since diagnosis',adj=0,cex.main=1)
dev.off()

# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
# Predefined follow-up to distinguish demographic composition shifts from era gaps.
# Uses all-period eligible CUP as a fixed measured-covariate reference, not a causal population.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(splines)})
setDTthreads(2)
out<-'outputs/stage2_preanalysis'
e<-readRDS('work/stage2_eligible.rds')
km<-function(z,w){q<-data.table(t=z$t,event=z$event,w=w)[,.(exit=sum(w),D=sum(w*event)),by=t][order(t)];q[,risk:=rev(cumsum(rev(exit)))];q[,S:=cumprod(pmax(0,1-D/risk))];s<-sapply(c(90,365.25),function(tau){j<-which(q$t<=tau);if(length(j))q$S[max(j)]else 1});qt<-q[t<365.25];c(OS90=s[1],OS365=s[2],RMST365=sum(diff(c(0,qt$t,365.25))*c(1,qt$S)))}
results<-list();diags<-list()
for(def in c('Broad','Strict8140')){
 target<-if(def=='Broad')e[group=='CUP']else e[group=='CUP'&hist==8140]
 target<-copy(target);target[,target_flag:=1L]
 for(er in levels(e$era))for(g in c('CUP','Lung','Pancreas','Colorectal','Stomach')){
  z<-copy(e[era==er&group==g&(def=='Broad'|hist==8140)]);z[,target_flag:=0L]
  zz<-rbind(target,z)
  terms<-c('ns(age,df=4)','old90','factor(sex)','factor(race)','factor(marital3)','factor(dc3)')
  terms<-terms[vapply(terms,function(v){if(grepl('^factor',v))uniqueN(zz[[sub('^factor\\((.*)\\)$','\\1',v)]])>1 else if(v=='old90')uniqueN(zz$old90)>1 else TRUE},logical(1))]
  fit<-glm(as.formula(paste('target_flag~',paste(terms,collapse='+'))),data=zz,family=binomial())
  stopifnot(fit$converged)
  p<-tail(fit$fitted.values,nrow(z));w<-p/(1-p)
  X<-model.matrix(fit)[,-1,drop=FALSE];nref<-nrow(target);Xt<-X[seq_len(nref),,drop=FALSE];Xz<-tail(X,nrow(z));den<-sqrt((apply(Xt,2,var)+apply(Xz,2,var))/2);smd<-(colMeans(Xt)-colSums(Xz*w)/sum(w))/den
  k<-km(z,w)
  results[[length(results)+1]]<-data.table(definition=def,era=er,group=g,n=nrow(z),target_n=nref,OS90=k[1],OS365=k[2],RMST365=k[3])
  diags[[length(diags)+1]]<-data.table(definition=def,era=er,group=g,n=nrow(z),ESS=sum(w)^2/sum(w^2),weight_max=max(w),max_abs_SMD=max(abs(smd),na.rm=TRUE))
 }
}
r<-rbindlist(results);dg<-rbindlist(diags)
fwrite(r,file.path(out,'17_fixed_target_era_OS.csv'),bom=TRUE)
fwrite(dg,file.path(out,'18_fixed_target_diagnostics.csv'),bom=TRUE)
contr<-merge(r[group!='CUP'],r[group=='CUP',.(definition,era,CUP_OS90=OS90,CUP_OS365=OS365,CUP_RMST365=RMST365)],by=c('definition','era'))
contr[,`:=`(OS90_diff_pp=100*(CUP_OS90-OS90),OS365_diff_pp=100*(CUP_OS365-OS365),RMST_diff_days=CUP_RMST365-RMST365)]
fwrite(contr,file.path(out,'19_fixed_target_era_contrasts.csv'),bom=TRUE)
png(file.path(out,'fixed_target_era_OS.png'),width=1600,height=900,res=150)
par(mar=c(6,5,4,12),xpd=FALSE)
cols<-c(CUP='#222222',Lung='#2684A6',Pancreas='#C57B19',Colorectal='#558B45',Stomach='#91539B')
plot(1:4,rep(0,4),type='n',ylim=c(0,65),xaxt='n',xlab='Diagnosis era',ylab='Standardized 1-year observed OS (%)',main='Fixed CUP covariate reference: preliminary point estimates')
axis(1,at=1:4,labels=levels(e$era));grid()
for(g in names(cols)){z<-r[definition=='Broad'&group==g];z[,i:=match(era,levels(e$era))];setorder(z,i);lines(z$i,100*z$OS365,type='b',col=cols[g],pch=19,lwd=2)}
par(xpd=NA)
legend(4.15,60,legend=names(cols),col=cols,lty=1,pch=19,bty='n',cex=.9)
mtext('Demographic/confirmation standardization only; no definitive confidence intervals or causal interpretation.',side=1,line=4.8,cex=.65)
dev.off()
print(r)
print(dg[,.(max_SMD=max(max_abs_SMD),min_ESS=min(ESS)),by=definition])

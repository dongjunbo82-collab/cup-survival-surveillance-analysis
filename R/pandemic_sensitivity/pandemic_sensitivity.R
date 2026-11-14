# Generated artifacts are written to local output directories.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/glandular_models', 'work/landmark',
            'outputs/standardization', 'outputs/cohort_checks',
            'outputs/definition_models', 'outputs/era_models', 'outputs/glandular_models',
            'outputs/pandemic_sensitivity', 'outputs/absolute_outcomes',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
library(data.table)
setDTthreads(1)
root <- Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'))
dest <- file.path(root,'outputs/pandemic_sensitivity');dir.create(dest,showWarnings=FALSE)
src <- file.path(root,'outputs/glandular_models')
pt <- fread(file.path(src,'04_point_estimates.csv'))
ss <- readRDS(file.path(root,'work/glandular_models/summary_arrays.rds'))
A<-ss$A;stopifnot(dim(A)[1]==2000,!anyNA(A))
ci<-function(v)c(n=length(v),lo=unname(quantile(v,.025)),hi=unname(quantile(v,.975)),lo_family=unname(quantile(v,.00625)),hi_family=unname(quantile(v,.99375)))
lab<-function(def,er,g)paste('Era',def,er,g,sep='__')
pval<-function(l)pt[label==l&measure=='OS365',estimate]
rows<-list()
for(def in c('Glandular','MinimalNoNE'))for(early in c('2004-2007','2008-2014'))for(g in c('Lung','Pancreas','Colorectal','Stomach')){
 ll<-c(lab(def,'2015-2019','CUP'),lab(def,'2015-2019',g),lab(def,early,'CUP'),lab(def,early,g))
 v<-A[,ll[1],'OS365']-A[,ll[2],'OS365']-A[,ll[3],'OS365']+A[,ll[4],'OS365']
 pp<-sum(vapply(ll,pval,numeric(1))*c(1,-1,-1,1))
 rows[[length(rows)+1]]<-cbind(data.table(definition=def,earlier=early,later='2015-2019',control=g,estimate=pp),as.data.table(as.list(ci(v))))
}
res<-rbindlist(rows);fwrite(res,file.path(dest,'33_prepandemic_gap_changes_CI.csv'),bom=TRUE)
# Check the first contrast against already preserved era contrasts.
reference_result<-fread(file.path(src,'17_gap_changes_CI.csv'))[era=='2015-2019'&measure=='OS365']
cmp<-merge(res[earlier=='2004-2007'],reference_result,by=c('definition','control'),suffixes=c('_comparison','_reference'))
stopifnot(nrow(cmp)==8,max(abs(cmp$estimate_comparison-cmp$estimate_reference))<1e-12,max(abs(cmp$lo_comparison-cmp$lo_reference))<1e-12,max(abs(cmp$hi_comparison-cmp$hi_reference))<1e-12)
dg<-fread(file.path(src,'05_weight_diagnostics.csv'))
annual<-pt[scope=='Annual'&definition=='Glandular'&measure=='OS365'&group%in%c('CUP','Lung','Pancreas')]
ac<-rbindlist(lapply(seq_len(nrow(annual)),function(i){z<-annual[i];cbind(z,as.data.table(as.list(ci(A[,z$label,'OS365']))))}))
ac<-merge(ac,dg[,.(label,support_pass,max_SMD,ESS)],by='label',all.x=TRUE)
stopifnot(nrow(ac)==57,all(ac$support_pass))
setorder(ac,group,period);fwrite(ac,file.path(dest,'34_annual_standardized_OS365_CI.csv'),bom=TRUE)
png(file.path(dest,'FigureS1_annual_survival.png'),width=2700,height=1650,res=300)
par(mfrow=c(1,3),mar=c(4.3,4.0,2.6,.8),oma=c(1.2,0,1,0),family='sans',cex=.95)
cols<-c(CUP='#202020',Lung='#2166AC',Pancreas='#D97706')
for(g in names(cols)){
 z<-ac[group==g][order(as.integer(period))];x<-as.integer(z$period)
 plot(NA,xlim=c(2004,2022),ylim=c(0,55),xlab='Diagnosis year',ylab='Standardized 1-year overall survival (%)',xaxt='n',yaxs='i')
 rect(2019.5,-1,2022.5,56,col='#F1F1F1',border=NA);abline(h=seq(0,50,10),col='#E6E6E6');abline(v=2019.5,lty=3,col='#999999')
 polygon(c(x,rev(x)),c(100*z$lo,rev(100*z$hi)),col=adjustcolor(cols[g],alpha.f=.16),border=NA)
 lines(x,100*z$estimate,col=cols[g],lwd=2);points(x,100*z$estimate,col=cols[g],pch=16,cex=.6)
 axis(1,at=c(2004,2008,2012,2016,2019,2022),labels=c('2004','2008','2012','2016','2019','2022'),las=2,cex.axis=.8)
 title(main=g,adj=0)
}
mtext('Bands: pointwise 95% bootstrap CI. Gray background: diagnosis years 2020-2022.',side=1,outer=TRUE,line=.2,cex=.8)
dev.off()
fwrite(data.table(check=c('2000 shared bootstrap draws','Paired contrast verified','57 annual models pass support'),passed=TRUE),file.path(dest,'35_sensitivity_QA.csv'))
writeLines(c('Exploratory contrasts use the shared bootstrap estimates.',
'Uses the preserved paired 2,000 patient bootstrap draws and fixed full-period pooled CUP target.',
'33: changes in CUP-minus-control 1-year OS gaps, 2015-2019 minus the named earlier era. Probability units; multiply by 100 for percentage points.',
'lo/hi: nominal 95% percentile CI. lo_family/hi_family: exploratory four-control Bonferroni intervals within each era contrast; not a separate prespecified primary family.',
'Both earlier-era comparisons use the same paired bootstrap draws.',
'Excluding pandemic-era diagnoses does not exclude all pandemic-era follow-up: some 2019 diagnoses have 1-year follow-up into 2020. This is not a causal estimate of the pandemic effect.',
'34: annual OS for CUP, lung and pancreas; all 57 fits pass support. The figure uses pointwise, not simultaneous, confidence bands.',
'Annual colorectal and stomach series omitted because their support criteria fail in some years.'),file.path(dest,'README.txt'))
print(res[definition=='Glandular'])

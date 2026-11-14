# Self-contained computational checks. Synthetic records are generated in memory.
suppressPackageStartupMessages({library(data.table);library(splines);library(survival)})
set.seed(1701)
n <- 800L
e <- data.table(age=sample(20:85,n,TRUE),old90=0L,sex=sample(c('Female','Male'),n,TRUE),
                race=sample(c('A','B'),n,TRUE),marital3=sample(c('Married','Other'),n,TRUE),
                dc3=sample(c('Histology','Cytology'),n,TRUE),
                group=rep(c('CUP','Lung'),each=n/2),era=rep('Early',n),
                t=sample(seq(.5,500,.5),n,TRUE),event=rbinom(n,1,.65))
e[,`:=`(days=t,months=t/(365.25/12))]
source('R/standardization/analysis_functions.R')
m <- build_model(which(e$group=='CUP'),which(e$group=='Lung'),'smoke','Synthetic','Early','Lung')
z <- fit_model(m,diagnostic=TRUE);stopifnot(identical(z$error,''),all(is.finite(z$values)))
# Compressed logistic likelihood versus an explicit stacked-record fit.
stack <- rbind(e[m$target],e[m$source]);stack[,role:=c(rep(1,length(m$target)),rep(0,length(m$source)))]
f <- glm(role~ns(age,df=4)+factor(sex)+factor(race)+factor(marital3)+factor(dc3),data=stack,family=binomial(),control=glm.control(maxit=40))
p <- tail(f$fitted.values,length(m$source));stopifnot(max(abs(z$odds-p/(1-p)))<1e-9)
# Weighted KM against survival::survfit, at the two reported horizons.
km <- survfit(Surv(e$t[m$source],e$event[m$source])~1,weights=z$odds)
ss <- summary(km,times=c(90,365.25),extend=TRUE)$surv
stopifnot(max(abs(ss-z$values[c('OS90','OS365')]))<1e-12)
stopifnot(abs(z$values['Conditional365given90']-ss[2]/ss[1])<1e-12)
# A common patient multiplicity is reused for source/reference roles.
counts <- tabulate(sample.int(n,n,TRUE),nbins=n)
zb <- fit_model(m,counts);stopifnot(identical(zb$error,''),all(is.finite(zb$values)))
unit <- fit_model(m,rep(1L,n));stopifnot(max(abs(unit$values-z$values))<1e-12)
cat('PASS: compressed likelihood, weighted KM, conditional ratio and shared multiplicities.\n')

# Generated artifacts are written to local output directories.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/glandular_models', 'work/landmark',
            'outputs/standardization', 'outputs/cohort_checks',
            'outputs/definition_models', 'outputs/era_models', 'outputs/glandular_models',
            'outputs/pandemic_sensitivity', 'outputs/absolute_outcomes',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(survival);library(splines)})
setDTthreads(2)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'));out<-file.path(root,'outputs/standardization/preparation_audit');dir.create(file.path(root,'work'),recursive=TRUE,showWarnings=FALSE)
src<-Sys.getenv('CUP_EXPORT_DIR');stopifnot(nzchar(src),dir.exists(src))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
save_table<-function(n,x)fwrite(x,file.path(out,paste0(n,'.csv')),bom=TRUE,na='NA')
files<-c(CUP='S2_CUP_Adeno_2004_2022_labels.txt',Lung='S2_Lung_Adeno_Distant_2004_2022_labels.txt',Pancreas='S2_Pancreas_Adeno_Distant_2004_2022_labels.txt',Colorectal='S2_Colorectal_Adeno_Distant_2004_2022_labels.txt',Stomach='S2_Stomach_Adeno_Distant_2004_2022_labels.txt')
newnames<-c('id','record','sequence','year','age_label','age_band','sex','race','marital','site','hist','behavior','dc','source','months','flag','days','vital','followyear','income','rural','stage','bone','brain','liver','lung')
inventory<-list();flows<-list();raws<-list();refnames<-NULL
for(g in names(files)){
 z<-fread(file.path(src,files[[g]]),colClasses='character',na.strings=NULL)
 if(is.null(refnames))refnames<-copy(names(z))
 stopifnot(ncol(z)==26,identical(names(z),refnames))
 setnames(z,newnames)
 for(v in names(z))set(z,j=v,value=trimws(z[[v]]))
 z[,group:=g]
 for(v in c('year','hist','months','days'))set(z,j=v,value=suppressWarnings(as.numeric(z[[v]])))
 z[months==9999,months:=NA_real_];z[days>=999999,days:=NA_real_]
 z[,age:=as.numeric(sub('^([0-9]+).*','\\1',age_label))]
 validhist<-with(z,(hist>=8140&hist<=8239)|(hist>=8250&hist<=8384)|(hist>=8440&hist<=8490)|hist%in%c(8550,8551,8570:8574,8576))
 siteok<-switch(g,CUP=grepl('^C80.9',z$site),Lung=grepl('^C34.',z$site),Pancreas=grepl('^C25.',z$site),Colorectal=grepl('^C18.|^C19.9|^C20.9',z$site),Stomach=grepl('^C16.',z$site))
 stopifnot(all(validhist),all(siteok),all(z$year%in%2004:2022),all(z$age>=20),all(z$behavior=='Malignant'))
 if(g!='CUP')stopifnot(all(z$stage=='Distant site(s)/node(s) involved'))
 inventory[[g]]<-data.table(group=g,filename=files[[g]],raw_n=nrow(z),fields=ncol(z)-2,year_min=min(z$year),year_max=max(z$year),duplicate_id_record=sum(duplicated(z[,.(id,record)])),md5=unname(tools::md5sum(file.path(src,files[[g]]))))
 mask<-rep(TRUE,nrow(z));fl<-data.table(group=g,step='Raw',remaining=nrow(z),excluded=0L)
 rules<-list('Sequence00/01'=z$sequence%in%c('One primary only','1st of 2 or more primaries'),'Exclude DCO/autopsy'=!tolower(z$source)%in%c('autopsy only','death certificate only'),'Microscopic1/2/4'=z$dc%in%c('Positive histology','Positive exfoliative cytology, no positive histology','Positive microscopic confirm, method not specified'),'Known survival and vital'=is.finite(z$days)&z$days>=0&is.finite(z$months)&z$months>=0&z$vital%in%c('Alive','Dead'))
 for(nm in names(rules)){nmask<-mask&rules[[nm]];fl<-rbind(fl,data.table(group=g,step=nm,remaining=sum(nmask),excluded=sum(mask)-sum(nmask)));mask<-nmask}
 z[,eligible:=mask];flows[[g]]<-fl;raws[[g]]<-z
}
save_table('01_source_audit',rbindlist(inventory))
save_table('02_selection_flow',rbindlist(flows))
d<-rbindlist(raws);rm(raws);gc()
save_table('03_raw_histology',d[,.N,by=.(group,hist)][order(group,-N)])
save_table('04_raw_stage',d[,.N,by=.(group,stage)])
e<-copy(d[eligible==TRUE]);e[,event:=as.integer(vital=='Dead')];e[,t:=pmax(days,0.5)]
e[,old90:=as.integer(age_band=='90+ years')]
e[,marital3:=fcase(grepl('^Married',marital),'Married',marital%in%c('Unknown',''),'Unknown',default='Not_married')]
e[,dc3:=fcase(dc=='Positive histology','Histology',grepl('cytology',dc),'Cytology',default='Micro_unspecified')]
e[,rural3:=fcase(grepl('^Counties in metropolitan',rural),'Metro',grepl('^Nonmetropolitan',rural),'Nonmetro',default='Unknown')]
e[,income_num:=suppressWarnings(as.numeric(gsub(',','',sub('.*?([0-9][0-9,]*).*','\\1',income))))]
e[grepl('^<',income),income_num:=0]
e[,income4:=fifelse(is.na(income_num),'Unknown',as.character(cut(income_num,c(-Inf,69999,89999,109999,Inf),labels=c('Below70k','70to90k','90to110k','110kplus'))))]
e[,era:=cut(year,c(2003,2007,2014,2019,2022),labels=c('2004-2007','2008-2014','2015-2019','2020-2022'))]
e[,organ_complete:=bone%in%c('Yes','No')&brain%in%c('Yes','No')&liver%in%c('Yes','No')&lung%in%c('Yes','No')]
e[,organ_positive:=bone=='Yes'|brain=='Yes'|liver=='Yes'|lung=='Yes']
saveRDS(e,file.path(root,'work/eligible_cohort.rds'))

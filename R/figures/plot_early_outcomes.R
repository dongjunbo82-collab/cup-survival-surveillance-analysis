# Generated artifacts are written to local output directories.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/glandular_models', 'work/landmark',
            'outputs/standardization', 'outputs/cohort_checks',
            'outputs/definition_models', 'outputs/era_models', 'outputs/glandular_models',
            'outputs/pandemic_sensitivity', 'outputs/absolute_outcomes',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(jsonlite)})
setDTthreads(1)
# Register the Arial-compatible metrics for the temporary base-PDF alignment device.
grDevices::pdfFonts(Arial=grDevices::pdfFonts('Helvetica')[[1]])
out<-'outputs/figures';work<-'work/landmark/figures';dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(work,recursive=TRUE,showWarnings=FALSE)
# Contract: a quantitative grid contrasts the early-death burden (A) with prognosis conditional on surviving 90 days (B); all estimates already exist.
x<-fread('outputs/early_landmark/four_era_absolute_probabilities.csv')[measure%in%c('Death90','Conditional365given90')]
stopifnot(nrow(x)==16,all(x$lo<=x$estimate&x$estimate<=x$hi),all(x$draws==2000))
x[,era:=factor(era,levels=c('2004-2007','2008-2014','2015-2019','2020-2022'))];x[,group:=factor(group,levels=c('CUP','Pancreas'))]
fwrite(x,file.path(out,'Figure2_source_data.csv'))
base<-theme_classic(base_size=8,base_family='Arial')+theme(plot.title=element_text(size=9,face='bold',hjust=0),axis.title=element_text(size=8),axis.text=element_text(size=7.5,colour='black'),axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=8),plot.margin=margin(6,7,6,6))
panel<-function(me,title,ylab){ggplot(x[measure==me],aes(x=era,y=estimate*100,colour=group,shape=group))+geom_errorbar(aes(ymin=lo*100,ymax=hi*100),width=.12,linewidth=.45,position=position_dodge(width=.3))+geom_point(size=2,position=position_dodge(width=.3))+scale_colour_manual(values=c(CUP='#202020',Pancreas='#D97706'),labels=c('CUP','Pancreatic cancer'))+scale_shape_manual(values=c(CUP=16,Pancreas=17),labels=c('CUP','Pancreatic cancer'))+scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20),expand=expansion(mult=c(0,.03)))+scale_x_discrete(labels=c('2004-2007','2008-2014','2015-2019','2020-2022'))+labs(title=title,x='Diagnosis era',y=ylab)+base}
p<-panel('Death90','A  Early mortality','90-day mortality (%)')+panel('Conditional365given90','B  Conditional survival','Survival to 1 year among 90-day survivors (%)')+plot_layout(ncol=2,guides='collect')&theme(legend.position='bottom')
w<-183/25.4;h<-105/25.4
ragg::agg_png(file.path(out,'Figure2_early_outcomes.png'),width=w,height=h,units='in',res=600);print(p);dev.off()
svglite::svglite(file.path(out,'Figure2_early_outcomes.svg'),width=w,height=h);print(p);dev.off()
ragg::agg_tiff(file.path(out,'Figure2_early_outcomes.tiff'),width=w,height=h,units='in',res=600,compression='lzw');print(p);dev.off()
ragg::agg_png(file.path(work,'Figure2_preview.png'),width=w,height=h,units='in',res=150);print(p);dev.off()
write_json(list(backend='R ggplot2/patchwork',width_mm=183,height_mm=105,source_rows=16,source_total_rows=24,selection='Two designated outcomes; OS365 remains in Table S20',interval='Existing pointwise nominal 95% bootstrap percentile intervals',new_models=FALSE,new_tests=FALSE,trend_lines=FALSE),file.path(work,'figure_contract.json'),pretty=TRUE)

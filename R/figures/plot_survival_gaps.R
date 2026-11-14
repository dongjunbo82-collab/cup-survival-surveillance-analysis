library(data.table);library(ggplot2);library(patchwork)
out<-'outputs/figures'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Claim: absolute survival changed against CUP, with the same gap-change directions
# under a narrower histologically confirmed target. Panel A establishes levels;
# panel B estimates paired gap changes and shows the restricted-target sensitivity.
aa<-fread('outputs/absolute_outcomes/36_absolute_OS365.csv')[era%in%c('2004-2007','2020-2022')]
aa[,group:=factor(group,levels=rev(c('CUP','Lung','Pancreas','Colorectal','Stomach')))];aa[,era:=factor(era,levels=c('2004-2007','2020-2022'))]
th<-theme_classic(base_size=8,base_family='sans')+theme(plot.title=element_text(size=9,face='bold',hjust=0),axis.title=element_text(size=8),legend.title=element_blank(),legend.position='bottom',legend.text=element_text(size=7),plot.margin=margin(7,7,7,7),axis.ticks.y=element_blank())
p1<-ggplot(aa,aes(x=estimate*100,y=group))+geom_line(aes(group=group),orientation='y',colour='grey65',linewidth=.5)+geom_errorbar(aes(xmin=lo*100,xmax=hi*100),orientation='y',width=.08,linewidth=.4)+geom_point(aes(shape=era),size=2.1)+scale_shape_manual(values=c(1,16))+scale_x_continuous(limits=c(0,65),breaks=seq(0,60,20))+labs(title='A  Absolute 1-year survival',x='Standardized 1-year survival (%)',y=NULL)+th
main<-fread('outputs/glandular_models/17_gap_changes_CI.csv')[definition=='Glandular'&era=='2020-2022'&measure=='OS365',.(control,estimate_pp=estimate*100,lo=lo_family*100,hi=hi_family*100,target='Main CUP reference')]
tissue<-fread('outputs/tissue_inference/gap_changes_CI.csv')[,.(control,estimate_pp,lo=lo_adjusted_pp,hi=hi_adjusted_pp,target='Histologically confirmed reference')]
bb<-rbind(main,tissue);bb[,control:=factor(control,levels=rev(c('Lung','Pancreas','Colorectal','Stomach')))];bb[,target:=factor(target,levels=c('Main CUP reference','Histologically confirmed reference'))]
pd<-position_dodge(width=.42)
p2<-ggplot(bb,aes(x=estimate_pp,y=control,colour=target,shape=target))+geom_vline(xintercept=0,linetype=2,colour='grey60',linewidth=.35)+geom_errorbar(aes(xmin=lo,xmax=hi),orientation='y',position=pd,width=.12,linewidth=.5)+geom_point(position=pd,size=2)+scale_colour_manual(values=c('#202020','#2166AC'))+scale_shape_manual(values=c(16,17))+scale_x_continuous(limits=c(-22,1),breaks=c(-20,-10,0))+labs(title='B  Change in CUP-minus-benchmark gap',x='1-year survival gap change (percentage points)',y=NULL)+th
fig<-p1+p2+plot_layout(widths=c(1,1),guides='keep')
ggsave(file.path(out,'Figure1_absolute_survival.png'),fig,width=8.2,height=4.3,units='in',dpi=450,bg='white',device=ragg::agg_png)
ggsave(file.path(out,'Figure1_absolute_survival.svg'),fig,width=8.2,height=4.3,units='in',device=svglite::svglite)
fwrite(bb,file.path(out,'Figure1_gap_source_data.csv'))

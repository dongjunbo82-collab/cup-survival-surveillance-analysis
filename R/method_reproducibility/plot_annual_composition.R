# Portable publication copy. Generated artifacts remain local and are ignored by Git.
root <- Sys.getenv('CUP_PROJECT_ROOT', unset=normalizePath(getwd(), winslash='/'))
setwd(root)
for (p in c('work', 'work/histology_review', 'work/cancer_control_revision',
            'outputs/method_reproducibility', 'outputs/stage2_preanalysis',
            'outputs/validation_stage3', 'outputs/final_audit', 'outputs/histology_revision',
            'outputs/pandemic_sensitivity', 'outputs/clinical_readability',
            'outputs/early_landmark', 'outputs/weight_capping', 'outputs/tissue_inference'))
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
library(data.table);library(jsonlite)
root<-Sys.getenv('CUP_PROJECT_ROOT',unset=normalizePath(getwd(),winslash='/'));out<-file.path(root,'outputs/method_reproducibility')
z<-fread(file.path(out,'50_annual_CUP_composition.csv'));stopifnot(identical(z$year,2004:2022),sum(z$n)==27715)
# Figure contract: descriptive composition changes, not evidence of reclassification.
# Quantitative grid; equal axes widths, separate percentage scales for common/unknown categories.
draw<-function(){
 par(mfrow=c(1,2),mar=c(4.4,4.6,3.1,1.1),oma=c(1.4,0,0,0),family='Arial',cex=.95,ps=12)
 bounds<-list()
 for(i in 1:2){
  vars<-if(i==1)c('histology_pct','strict8140_pct')else c('marital_unknown_pct','race_unknown_pct')
  labs<-if(i==1)c('Histologic confirmation','Adenocarcinoma NOS 8140/3')else c('Unknown marital status','Unknown race/ethnicity')
  plot(NA,xlim=c(2004,2022),ylim=if(i==1)c(60,100)else c(0,10),xaxs='i',yaxs='i',xaxt='n',xlab='Diagnosis year',ylab='CUP records (%)')
  abline(h=if(i==1)seq(60,100,10)else seq(0,10,2),col='#EEEEEE');axis(1,at=c(2004,2008,2012,2016,2020,2022),labels=c(2004,2008,2012,2016,2020,2022),las=2,cex.axis=.85)
  for(j in 1:2){lines(z$year,z[[vars[j]]],col=c('#202020','#2166AC')[j],lty=j,lwd=2);points(z$year,z[[vars[j]]],pch=c(16,17)[j],col=c('#202020','#2166AC')[j],cex=.6)}
  title(main=if(i==1)'A  Diagnostic confirmation and morphology'else'B  Unknown recorded characteristics',adj=0,cex.main=.9)
  legend(if(i==1)'bottomleft'else'topleft',legend=labs,col=c('#202020','#2166AC'),lty=1:2,pch=c(16,17),bty='n',cex=.8)
  f<-par('fig');p<-par('plt');ds<-dev.size('in')*72
  bounds[[i]]<-list(id=c('a','b')[i],bbox_pt=c((f[1]+p[1]*diff(f[1:2]))*ds[1],(f[3]+p[3]*diff(f[3:4]))*ds[2],(f[1]+p[2]*diff(f[1:2]))*ds[1],(f[3]+p[4]*diff(f[3:4]))*ds[2]))
 }
 mtext('Unweighted annual proportions in the reviewed CUP cohort. Lines join observed years; no smoothing or trend test.',side=1,outer=TRUE,line=.25,cex=.78)
 list(schema_version=1,figure=list(width_pt=dev.size('in')[1]*72,height_pt=dev.size('in')[2]*72),panels=bounds,row_groups=list(list(id='composition',panels=c('a','b'))))
}
png(file.path(out,'FigureS3_annual_CUP_composition.png'),width=2400,height=1440,res=300);geo<-draw();dev.off()
svglite::svglite(file.path(out,'FigureS3_annual_CUP_composition.svg'),width=8,height=4.8);draw();dev.off()
write_json(geo,file.path(out,'FigureS3.alignment-layout.json'),auto_unbox=TRUE,pretty=TRUE)
# Original pooled-reference knot positions for every era model, held fixed in bootstrap.
e<-readRDS(file.path(root,'work/stage2_eligible.rds'));cfg<-fromJSON(file.path(root,'config/histology_definitions.json'));e<-e[hist%in%cfg$reviewed_keep_codes];ref<-e[group=='CUP']
knots<-list();for(er in levels(e$era))for(g in c('CUP','Lung','Pancreas','Colorectal','Stomach')){a<-splines::ns(c(ref$age,e[group==g&era==er,age]),df=4);knots[[length(knots)+1]]<-data.table(era=er,group=g,internal_knots=paste(attr(a,'knots'),collapse=';'),boundary_knots=paste(attr(a,'Boundary.knots'),collapse=';'))}
fwrite(rbindlist(knots),file.path(out,'53_primary_spline_settings.csv'))

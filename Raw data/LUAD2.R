if (T) {
 
  dir.create("Files")
  dir.create("data")
  dir.create("Figures")
  dir.create("00_origin_datas/GEO",recursive = T)
  dir.create("00_origin_datas/TCGA")
  dir.create("00_pre_datas/GEO",recursive = T)
  dir.create("00_pre_datas/TCGA")
}
library(cols4all)
library(clusterProfiler)
library(ggpubr)
library(stringr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
#library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
#library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)

get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=genes,model=mult_results))
}
my_mutiviolin=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                       #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                       bw=T,xlab='',ylab='score',title='',size=3,angle = 45, hjust = 1,
                       legend.position='top',fill='group',notch=F){
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(dat.melt,aes(x=type, y=value,fill=Group)) +
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill='Group')+
    geom_violin(trim = F,position=position_dodge(0.8))+  
    scale_fill_manual(values = group_cols)+
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust))
  return(p)
}
plotMutiBar <-plotMutiBar <-function(dat=Age1_compare,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T,color){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) 
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_npg()+theme(legend.position = "bottom")
  pg=pg+ggsci::scale_fill_npg()+scale_fill_manual(values = color, 
                                                  breaks = paste0('R', 1:nrow(dat)), 
                                                  labels = lbr, 
                                                  name = legTitle) 
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") +xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(CHl-Squared p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

wb_beeswarm_plot <- function(dat = NULL,
                             show_compare = T,
                             xlab = 'Groups',
                             ylab = '',
                             method = c('t.test', 'wilcox.test')[1],
                             col = mycolor,
                             leg.pos = c('top','left','right','bottom','none')[1],
                             title = NULL,
                             group = 'Cluster') {
  library(ggbeeswarm)
  colnames(dat) <- c('Feature', 'Cluster')
  
  
  p1 <- ggplot(dat, aes(Cluster, Feature, color = Cluster)) + geom_quasirandom(method = "frowney") +
    ggtitle(title) + scale_color_manual(values = col[1:length(unique(dat$Cluster))]) +
    xlab(xlab) + ylab(ylab) + guides(color=guide_legend(title = group)) + theme_classic() +
    theme(legend.position=leg.pos)
  
  
  if(show_compare){
    uni.group = as.character(unique(dat$Cluster))
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,
                                     method = method,
                                     label= "p.signif", 
                                     step_increase = 0.0)
  }
  return(p1)
}
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin(trim = F)+  
    scale_fill_manual(values = group_cols)+
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}
my_boxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                    fill= "Group",label=c("p.format",'p.signif')[1],notch=F,
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=ARGs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +
    scale_fill_manual(values = group_cols)+   #
    # if(length(names(table(group)))>2){
    #   test_method=''
    # }
    ggpubr::stat_compare_means(aes(group=Group), label = label, method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 
  return(p)
}

my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        test_method=c('t.test','wilcox.test','anova','kruskal.test')[2],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    # theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","grey"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(
      legend.title = element_blank(),#
      legend.position = leg.pos,
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#padj<0.05
  return(p)
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}

coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
wb_boxplot <- function(dat = data,
                       groups = groups,
                       xlab = '',
                       ylab = '',
                       xangle = 90,
                       title = 'Groups',
                       col = mycolor) {
  tmp.dat <- data.frame()
  for (ge in rownames(dat)) {
    print(ge)
    tmp <- data.frame(Samples = colnames(dat),
                      Genes = ge,
                      Values = as.numeric(dat[ge, ]),
                      Groups = groups)
    tmp.dat <- rbind(tmp.dat, tmp)
  }
  
  
  print(head(tmp.dat))
  library(ggpubr)
  if (length(unique(groups)) > 2) {
    tmp_plot <- ggplot(tmp.dat, 
                       aes(x=Genes, y=Values, fill=Groups)) +
      geom_boxplot(notch = F) +  
      stat_compare_means(method = "anova", label = "p.signif") +
      scale_fill_manual(values = col) + theme_classic() +
      theme(axis.text.x = element_text(angle=xangle, 
                                       hjust = 0.95,
                                       vjust = 0.95)) +
      xlab(xlab) + ylab(ylab) + guides(fill = guide_legend(title = title))
    return(tmp_plot)
  } else {
    tmp_plot <- ggplot(tmp.dat, 
                       aes(x=Genes, y=Values, fill=Groups)) +
      geom_boxplot(notch = F) +  
      stat_compare_means(method = "t.test", label = "p.signif") +
      scale_fill_manual(values = col) + theme_classic() +
      theme(axis.text.x = element_text(angle=xangle, 
                                       hjust = 0.95,
                                       vjust = 0.95)) +
      xlab(xlab) + ylab(ylab) + guides(fill = guide_legend(title = title))
    return(tmp_plot)
  }
}

####TCGA-LUAD################
tcga.cli1<-read.delim('00_origin_datas/TCGA/Merge_LUAD_clinical.txt',sep='\t',header = T)
colnames(tcga.cli1)[1:20]
tcga.cli1$number_pack_years_smoked
table(tcga.cli1$tobacco_smoking_history)
tcga.cli1=data.frame(Samples=tcga.cli1$A0_Samples,
                     Age=tcga.cli1$A17_Age,
                     Gender=tcga.cli1$A18_Sex,
                     Smoking_history=tcga.cli1$tobacco_smoking_history,
                     T.stage=tcga.cli1$A3_T,
                     N.stage=tcga.cli1$A4_N,
                     M.stage=tcga.cli1$A5_M,
                     Stage=tcga.cli1$A6_Stage)
tcga.cli1$Samples=paste0(tcga.cli1$Samples,'-01')
rownames(tcga.cli1)=tcga.cli1$Samples
head(tcga.cli1)
table(tcga.cli1$Smoking_history)

table(tcga.cli1$T.stage)
tcga.cli1$T.stage=gsub('[ab]','',tcga.cli1$T.stage)
tcga.cli1$T.stage[tcga.cli1$T.stage=='TX']<-NA

table(tcga.cli1$N.stage)
tcga.cli1$N.stage[tcga.cli1$N.stage=='NX'|tcga.cli1$N.stage=='']<-NA

table(tcga.cli1$M.stage)
tcga.cli1$M.stage=gsub('[ab]','',tcga.cli1$M.stage)
tcga.cli1$M.stage[tcga.cli1$M.stage=='MX'|tcga.cli1$M.stage=='']<-NA

table(tcga.cli1$Stage)
tcga.cli1$Stage=gsub('[AB]','',tcga.cli1$Stage)
tcga.cli1$Stage[tcga.cli1$Stage=='']<-NA
tcga.cli1$Stage=gsub('Stage ','',tcga.cli1$Stage)


tcga.pancancer.cli=read.xlsx
head(tcga.pancancer.cli)
tcga.cli2=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='LUAD'),]
head(tcga.cli2)
tcga.cli2=data.frame(Samples=paste0(tcga.cli2$bcr_patient_barcode,'-01'),
                     tcga.cli2[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
head(tcga.cli2)
tcga.cli2$OS.time
tcga.cli2=tcga.cli2 %>% drop_na(OS.time)
####
tcga.cli2=tcga.cli2[tcga.cli2$OS.time>0,]
dim(tcga.cli2)

tcga.cli=merge(tcga.cli1,tcga.cli2,by='Samples')
rownames(tcga.cli)=tcga.cli$Samples
tcga.cli=as.data.frame(tcga.cli)
fivenum(as.numeric(tcga.cli$Age))
tcga.cli$Age1=ifelse(as.numeric(tcga.cli$Age)>66,'>66','<=66')
dim(tcga.cli)
# 509  17
head(tcga.cli)


#
tcga_data<-read.delim('00_origin_datas/TCGA/Merge_RNA_seq_FPKM _LUAD.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga_data[1:4,1:4]
table(substr(colnames(tcga_data),14,15))
tcga_data <- exp_ensg2symbol(tcga_data)

sample_T=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==1)]#
sample_N=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==11)]#
length(sample_N)
length(sample_T)
tcga_type=data.frame(Samples=c(sample_T,sample_N),type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$type)
# Normal  Tumor 
# 59    513

genecode=read.delim('data/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]


range(tcga_data)
tcga.exp.all=log2(tcga_data[intersect(rownames(tcga_data),mrna_genecode$SYMBOL),tcga_type$Samples]+1)
range(tcga.exp.all)
tcga.exp=tcga.exp.all[,intersect(tcga.cli$Samples,sample_T)]
dim(tcga.exp)
# 19503   500
tcga.cli=tcga.cli[intersect(tcga.cli$Samples,sample_T),]
dim(tcga.cli)
saveRDS(tcga.exp.all,file = '00_pre_datas/TCGA/LUAD_FPKM_tcga.exp.all.RDS')
saveRDS(tcga.exp,file = '00_pre_datas/TCGA/LUAD_FPKM_tcga.exp.RDS')
########GSE31210#############
GSE31210 <- getGEOExpData('GSE31210')
saveRDS(GSE31210,file = "00_origin_datas/GEO/GSE31210.RDS")
#load('00_origin_datas/GEO/GSE31210.RData')

GSE31210.cli=GSE31210$Sample
GSE31210.cli=data.frame(Samples=GSE31210.cli$Acc,
                        Age=GSE31210.cli$`age (years)`,
                        Gender=GSE31210.cli$gender,
                        Status=GSE31210.cli$death,
                        OS.time=GSE31210.cli$`days before death/censor`)
rownames(GSE31210.cli)=GSE31210.cli$Samples
table(GSE31210.cli$Status)
GSE31210.cli=GSE31210.cli[which(GSE31210.cli$Status!='NULL'),]
GSE31210.cli$OS=ifelse(GSE31210.cli$Status=='alive',0,1)
GSE31210.cli <- GSE31210.cli[GSE31210.cli$OS.time>0,]
range(GSE31210.cli$OS.time)


GSE31210.exp=GSE31210$Exp$GPL570_54675_Data_col1
range(GSE31210.exp)
GSE31210.exp=log2(GSE31210.exp+1)
GSE31210.exp=exp_probe2symbol_v2(GSE31210.exp,GPL ='GPL570' )
range(GSE31210.exp)
dim(GSE31210.exp)
dim(GSE31210.cli)
GSE31210.exp=GSE31210.exp[,GSE31210.cli$Samples]
dim(GSE31210.exp)
#20549   226


####01.##########
dir.create('01_GBP_gene')
GBP.gene <- read.table("01_GBP_gene/GBP_Gene.txt",sep = "\t")
GBP.gene <- GBP.gene$V1
length(GBP.gene)
#7
GBP.gene <- intersect(GBP.gene,rownames(tcga.exp))

##################1.3############
GBP.score=ssGSEAScore_by_genes(tcga.exp.all,genes = GBP.gene)
GBP.score <- as.data.frame(t(GBP.score))
p1a <- my_violin(GBP.score ,group=tcga_type$type,ylab='GBP scores')
ggsave('01_GBP_gene/p1a.pdf',height = 6,width = 6)

####02WGCNA###########
dir.create("02_WGCNA")

library(WGCNA)
allowWGCNAThreads(nThreads = 36)#
enableWGCNAThreads(nThreads = 36)# 

my_mad <- function(x){mad(x,na.rm = TRUE)} #
wgcna_exp=t(tcga.exp)
m.mad <- apply(wgcna_exp,2,my_mad)
dim(tcga.exp)
# 18448   500
#
tpm_T2 <- wgcna_exp[,which(m.mad >max( quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]

#tpm_T2=tcga.exp.HBV[which(apply(tcga.exp.HBV,1,sd)>0.5),]
#tpm_T2=(2^tpm_T2-1)
range(tpm_T2)
pdf('02_WGCNA//wgcna.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()

tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.2,
                                 minModuleSize=60)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))
#18
write.csv(tpm_T2.module$Modules,file = "02_WGCNA//WGCNA_Modules.csv",col.names = T,row.names = T)
pdf('02_WGCNA/2.pdf',height = 6,width = 10)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
writeMatrix(tpm_T2.module$Modules,outpath = '02_WGCNA/tcga.wgcna.module.genes.txt')

pdf('02_WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()

#### 
# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('02_WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module genes",xlab = "", sub = "")
dev.off()


############### 
tcga_cli_use <-data.frame(tcga.cli[tcga.cli$Samples,c(1,8)],GBP.score=GBP.score[tcga.cli$Samples,])
head(tcga_cli_use)
tcga_cli_use=as.data.frame(tcga_cli_use[,-c(1,2)])
colnames(tcga_cli_use)="GBP.score"
tcga_cli_use.part=tcga_cli_use
str(tcga_cli_use.part)
#tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))


spms=tcga_cli_use.part
MEs_col<-tpm_T2.module$MEs
dim(MEs_col)
modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms[]
                  ,use = 'pairwise.complete.obs')

modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])
textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('02_WGCNA/5.pdf',width = 6,height =12)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#
geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))
#
geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))
head(geneTraitSignificance)
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames


module="greenyellow"
column= match(module, modNames)


moduleGenes <- c(tpm_T2.module$Modules[,'mergedColors']==module)
greenyellow.gene=names(which(moduleGenes))
length(greenyellow.gene)
# 422
writeMatrix(greenyellow.gene,'02_WGCNA/greenyellow.gene.txt',header = F)
table(tpm_T2.module$Modules[,'mergedColors'])

pdf('02_WGCNA/6.pdf',height = 6,width = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, ]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for immunity",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = 'red',lwd=2)
dev.off()



####03_###########
dir.create('03_diff_genes')
tcga.limma=mg_limma_DEG(exp =tcga.exp.all[,tcga_type$Samples],
                        group=tcga_type$type,ulab='Tumor',dlab = 'Normal')
tcga.limma$Summary

#######
tcga.Normal_Tumor.degs=tcga.limma$DEG
tcga.Normal_Tumor.degs=tcga.Normal_Tumor.degs[abs(tcga.Normal_Tumor.degs$logFC)>log2(1.5) & tcga.Normal_Tumor.degs$adj.P.Val<0.05,]
dim(tcga.Normal_Tumor.degs)
write.csv(tcga.Normal_Tumor.degs,'03_diff_genes/tcga.Normal_Tumor.degs.csv')

######3.1####
p_cutoff <-0.05 
fc_cutoff <- log2(1.5)
degs_dat=tcga.limma$DEG
degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                            ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))

#library(ggbreak)
library(ggplot2)
library(ggprism)

col=c("hotpink","sienna","grey")
ylab='-log10 (adj.PVal)'
xlab='log2 (FoldChange)'
leg.pos='right'
plot2a<- ggplot(degs_dat, aes(x=logFC, y=-log10(adj.P.Val), color=type)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=col) +
  theme_bw() +
  theme(legend.position = leg.pos) +
  ylab(ylab) +
  xlab(xlab) +
  geom_vline(xintercept=c(-fc_cutoff,fc_cutoff), lty=3, col="black", lwd=0.5) +
  geom_hline(yintercept = -log10(p_cutoff), lty=3, col="black", lwd=0.5) 
  #coord_cartesian(ylim=c(0, 25)) # 
  #scale_y_break(c(25,100),#
                # space = 0.3,#
                # scales = 0.5)+#
  #theme_prism(palette = "black_and_white",
              # base_fontface = "plain", 
              # base_family = "serif", 
              # base_size = 16,
              # base_line_size = 0.8,
              # axis_text_angle = 0)
ggsave('03_diff_genes/fig3a.pdf',plot2a,height = 6,width = 6)
######3.2#######
venny.gene=intersect(greenyellow.gene,rownames(tcga.Normal_Tumor.degs))
vennydata=list(greenyellow_module=greenyellow.gene,TCGA_DEG=rownames(tcga.Normal_Tumor.degs))
library("VennDiagram")
pdf('03_diff_genes/venn.pdf',height = 7,width = 6)
# 
venn.plot <- venn.diagram(
  x = vennydata,
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "plum2"),  # 
  alpha = 0.90,  # 
  cat.cex = 1.2,  # 
  cat.col = "gray10",  # 
  cat.fontface = "bold"  # 
)

# 
grid.draw(venn.plot)
dev.off()
######3.2 ########
length(venny.gene)
venny.gene.enrichment=mg_clusterProfiler(genes = venny.gene)
head(venny.gene.enrichment$GO_BP)
p3=list()
p3[[1]]=barplot(venny.gene.enrichment$GO_BP)
p3[[2]]=barplot(venny.gene.enrichment$KEGG)
p3cd=mg_merge_plot(p3,labels = c('C','D'),ncol = 2,widths = c(1.8,1))


write.xlsx(
           list(GO_BP=venny.gene.enrichment$GO_BP,
                KEGG=venny.gene.enrichment$KEGG),'03_diff_genes/venny.gene.enrichment.xlsx',
                overwrite = T)

savePDF('03_diff_genes/Fig3CD.pdf',p3cd,height = 10,width = 25)

####04.########
dir.create('04_model')
##########

sig.gene.cox=cox_batch(dat = tcga.exp[intersect(venny.gene,rownames(tcga.exp)),tcga.cli$Samples],
                        time = tcga.cli$OS.time,event = tcga.cli$OS)
sig.gene.cox



table(sig.gene.cox$p.value<0.01)
# FALSE  TRUE 
#   71    21 
pre.genes=rownames(sig.gene.cox[sig.gene.cox$p.value<0.01,])
length(pre.genes)#21
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
write.csv(sig.gene.cox,'04_model/sig.cox.csv')

##########
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(pre.genes,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
"-0.256*CBFA2T3+0.469*CD19+-0.381*MS4A1+-0.464*P2RX1"
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)

###########
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[names(lan), tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
fig_ggforest <- survminer::ggforest(cox,data=tcga_model_data)
ggsave("04_model/ggforest.pdf",height =4,width =7 )

####
risktype.col=c('#FF5154',"#704214")

risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.subtype.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=risk.tcga)
######KM####
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(tcga.risktype.cli$Riskscore),'High','Low')

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,3,5))
tcga.roc
tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                  data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = T,pval = T,risk.table = T, 
                      fun = "pct",size = 1,surv.median.line = 'hv',
                      title='TCGA-LUAD',legend.title='Risktype',
                      legend.labs = c('High','Low'),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,
                      ylab='Overall Survival(OS)',
                      legend=c(0.85,0.8),#
                      ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(3,1),align = 'v')

tcga.km.OS

tcga.risktype.cli$Status=ifelse(tcga.risktype.cli$OS==0,'Alive','Dead')
tcga.model.p=my_riskplot(cli_dat = tcga.risktype.cli,cols =risktype.col,xlab = 'sample',
            a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = median(tcga.risktype.cli$Riskscore),labs = '')


tcga.km.DSS=ggsurvplot(fit=survfit(Surv(DSS.time/365, DSS) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',
                       title='TCGA-LUAD',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,
                       ylab='Disease-Specific Survival(DSS)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.DSS=mg_merge_plot(tcga.km.DSS$plot,tcga.km.DSS$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.DSS

tcga.km.DFI=ggsurvplot(fit=survfit(Surv(DFI.time/365, DFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',
                       title='TCGA-LUAD',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,
                       ylab='Disease-Free interval(DFI)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.DFI=mg_merge_plot(tcga.km.DFI$plot,tcga.km.DFI$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.DFI

tcga.km.PFI=ggsurvplot(fit=survfit(Surv(PFI.time/365, PFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',title='TCGA-LUAD',
                       legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       palette = risktype.col,
                       ylab='Progression-Free Interval(PFI)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.PFI=mg_merge_plot(tcga.km.PFI$plot,tcga.km.PFI$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.PFI
###########
tcga.barplot=my_mutibarplot(df=table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')
########
model.gene.df=data.frame(tcga.cli[,c('OS','OS.time')],t(tcga.exp[names(lan),tcga.cli$Samples]))
head(model.gene.df)
module.gene.km=list()
for (i in 1:length(names(lan))) {
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km[[i]]=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                             data = model.gene.df1),
                                 data=model.gene.df1,
                                 conf.int = F,pval = T,risk.table = T, 
                                 fun = "pct",size = 1,surv.median.line = 'hv',
                                 title='TCGA',legend.title=names(lan)[i],
                                 # legend.labs = c('High','Low'),
                                 linetype = c("solid", "dashed","strata")[1],
                                 palette = risktype.col,
                                 ylab='Overall Survival(OS)',
                                 legend=c(0.8,0.8),#
                                 ggtheme = theme_bw(base_size = 12))
  module.gene.km[[i]]=module.gene.km[[i]]$plot
}
tcga.module.km=mg_merge_plot(module.gene.km,ncol=3,nrow=2)



########4.2 ##########
GSE31210_model_data <- data.frame(GSE31210.cli[, c("OS.time", "OS")],
                             t(GSE31210.exp[intersect(names(lan),rownames(GSE31210.exp)), GSE31210.cli$Samples]))
colnames(GSE31210_model_data) <- gsub('-', '_', colnames(GSE31210_model_data))

risk.GSE31210=as.numeric(lan%*%as.matrix(t(GSE31210_model_data[GSE31210.cli$Samples,names(lan)])))

GSE31210.risktype.cli=data.frame(GSE31210.cli,Riskscore=risk.GSE31210)
GSE31210.risktype.cli$Risktype=ifelse(GSE31210.risktype.cli$Riskscore>median(risk.GSE31210),'High','Low')
GSE31210.roc=ggplotTimeROC(GSE31210.risktype.cli$OS.time,
                           GSE31210.risktype.cli$OS,
                           GSE31210.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
GSE31210.roc
GSE31210.km=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                   data = GSE31210.risktype.cli),
                       data=GSE31210.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='GSE31210',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ylab='Overall Survival(OS)',
                       legend=c(0.85,0.25),#
                       ggtheme = theme_bw(base_size = 12))
GSE31210.km=mg_merge_plot(GSE31210.km$plot,GSE31210.km$table,nrow=2,heights = c(3,1),align = 'v')
GSE31210.km
#####GSE31210###
gse31210.model.p=my_riskplot(cli_dat = GSE31210.risktype.cli,cols = risktype.col,xlab = 'sample',
            a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = median(GSE31210.risktype.cli$Riskscore),labs = '')

########
my_mutibarplot=function(df,xlab='group',leg.title='',cols=pal_d3()(10)[5:6]){
  prop.pval=round(chisq.test(df)$p.value,2)#round(-log10(chisq.test(df)$p.value),2)
  if( prop.pval<0.001)
    prop.pval='<0.001'
  df.prop=prop.table(df,margin=2)
  df.prop=reshape2::melt(df.prop)
  colnames(df.prop)<-c("type","group","Percentage")
  df.prop$Percentage<-round(df.prop$Percentage,digits=2)
  p=ggplot(df.prop,aes(x=group,y=Percentage,fill=type))+
    geom_bar(position = "fill",stat="identity")+
    scale_fill_manual(values = cols)+
    xlab(xlab)+labs(fill = leg.title,title = 'Chi-Squared Test',subtitle  =  paste0('pvalue  ',prop.pval))+
    theme_bw()+theme(text=element_text(family = 'Times'),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
  p
  return(p)
}
GSE31210.barplot=my_mutibarplot(df=table(GSE31210.risktype.cli$Status,GSE31210.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')

##########
model.gene.df=data.frame(GSE31210.cli[,c('OS','OS.time')],t(GSE31210.exp[names(lan),GSE31210.cli$Samples]))
head(model.gene.df)
module.gene.km=list()
for (i in 1:length(names(lan))) {
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km[[i]]=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                             data = model.gene.df1),
                                 data=model.gene.df1,
                                 conf.int = F,pval = T,risk.table = T, 
                                 fun = "pct",size = 1,surv.median.line = 'hv',
                                 title='GSE31210',legend.title=names(lan)[i],
                                 # legend.labs = c('High','Low'),
                                 linetype = c("solid", "dashed","strata")[1],
                                 palette = risktype.col,
                                 ylab='Overall Survival(OS)',
                                 legend=c(0.7,0.35),#
                                 ggtheme = theme_bw(base_size = 12))
  module.gene.km[[i]]=mg_merge_plot(module.gene.km[[i]]$plot,module.gene.km[[i]]$table,nrow=2,heights = c(3,1),align = 'v')
}
GSE31210.module.km=mg_merge_plot(module.gene.km,ncol=3,nrow=2)
ggsave("04_model/GSE31210.module.km.pdf",height =10 ,width = 10)
###########
fig4_1=mg_merge_plot(mg_merge_plot(fig_ggforest,tcga.roc,labels = c('A','B'),widths = c(1.8,1)),
                   mg_merge_plot(tcga.km.DSS,tcga.km.PFI,tcga.km.OS,ncol=3,labels = LETTERS[3:5]),
                   mg_merge_plot(tcga.barplot,GSE31210.roc,GSE31210.km,ncol=3,labels = LETTERS[6:8]),
                   nrow=3)


ggsave('04_model/Fig4.pdf',fig4_1,height = 20,width = 20)



####05.##########
dir.create('05_Risktype.immune')
#  #############

#######estimate####
tcga_estimate <- immu_estimate(exp = tcga.exp)

p5a <-  my_mutiboxplot(tcga_estimate ,group = tcga.risktype.cli$Risktype,group_cols=risktype.col,legend.pos='top',ylab = "Estimate score")



#####TIMER#####
tcga.timer=immu_timer(tcga.exp)
p5b<- mg_PlotMutiBoxplot(data = tcga.timer[tcga.risktype.cli$Samples,],group =tcga.risktype.cli$Risktype,
                         legend.pos = 'top',group_cols = risktype.col,add = 'boxplot',test_method = 'wilcox.test',ylab = 'Timer score')

######ssgsea####
taga.ssgsea=immu.ssgsea(exp =tcga.exp,isTCGA = T )
saveRDS(taga.ssgsea,file ='05_Risktype.immune/tcga.immu.ssgsea.RDS')
tcga.immu.ssgsea <- readRDS("05_Risktype.immune/tcga.immu.ssgsea.RDS")
p5c <- groupViolin(tcga.immu.ssgsea [tcga.risktype.cli$Samples, ],
                                  tcga.risktype.cli$Risktype,
                                      ylab = 'ssgsea Immune Score',
                                      # legend.pos = 'br',
                                      group_col=risktype.col)




# ######
p5ab<- cowplot::plot_grid(p5a,
                          p5b,
                          ncol = 2,labels = c("A",'B'),rel_widths = c(1,1.8),align ="h" )

ggsave('05_Risktype.immune/p5ab.pdf',p5ab,height = 6,width = 15)

p5abc <- cowplot::plot_grid(p5ab,
                            p5c,
                          nrow = 2,labels = c("",'C'))

ggsave('05_Risktype.immune/p5ABC.pdf',p5abc,height = 12,width = 15)



####06#########
dir.create('06_nomo')
tcga.risktype.cli.use=tcga.risktype.cli[order(tcga.risktype.cli$Riskscore),]

colnames(tcga.risktype.cli.use)
cli.colors=list()
for(i in c(6:9,18,20)){
  var=tcga.risktype.cli.use[i]
  var.clas=names(table(tcga.risktype.cli.use[,i]))
  var.color=c4a('tableau.classic10',9)[1:length(var.clas)]
  names(var.color)=var.clas
  cli.colors=c(cli.colors,list(var.color))
}
names(cli.colors)=colnames(tcga.risktype.cli.use)[c(6:9,18,20)]


cli.colors[[6]]=risktype.col
names(cli.colors[[6]])=c('High','Low')
library(ComplexHeatmap)

tcga.riskscore.barplot=columnAnnotation(Riskscore = anno_barplot(tcga.risktype.cli.use$Riskscore,
                                                                 baseline =median(risk.tcga),
                                                                 bar_width=2,
                                                                 gp=gpar(fill=ifelse(tcga.risktype.cli.use$Riskscore>median(risk.tcga),
                                                                                     risktype.col[1],risktype.col[2])))
                                        ,annotation_name_side ='left'
                                        ,height =unit(4,'inches'))
draw(tcga.riskscore.barplot)

tcga.cli.heatmap=columnAnnotation(df = tcga.risktype.cli.use[,c(6:9,18,20)]
                                  ,annotation_name_side='left'
                                  ,annotation_height =unit(4,'inches')
                                  ,col = cli.colors)

ht_list=tcga.riskscore.barplot %v% tcga.cli.heatmap
pdf('06_nomo/Fig6a.pdf',height = 8,width = 12,onefile = F)
ht_list
dev.off()
###########
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
# table(tcga_cox_datas$Age1)
table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'

table(tcga_cox_datas$N.stage)
tcga_cox_datas$N.stage[tcga_cox_datas$N.stage=='N1'|tcga_cox_datas$N.stage=='N2'|tcga_cox_datas$N.stage=='N3']<-'N1+N2+N3'

table(tcga_cox_datas$M.stage)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'
#Age
colnames(tcga_cox_datas)
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

#T.stage
T.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.stage,
                                 data=tcga_cox_datas))
T.stage_sig_cox_dat <- data.frame(Names=rownames(T.stage_sig_cox[[8]]),
                                  HR = round(T.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(T.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(T.stage_sig_cox[[8]][,4],3),
                                  p.value=round(T.stage_sig_cox[[7]][,5],3))
T.stage_sig_cox_dat

#N.stage
N.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.stage,
                                 data=tcga_cox_datas))
N.stage_sig_cox_dat <- data.frame(Names=rownames(N.stage_sig_cox[[8]]),
                                  HR = round(N.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(N.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(N.stage_sig_cox[[8]][,4],3),
                                  p.value=round(N.stage_sig_cox[[7]][,5],3))
N.stage_sig_cox_dat

#M.stage
M.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.stage,
                                 data=tcga_cox_datas))
M.stage_sig_cox_dat <- data.frame(Names=rownames(M.stage_sig_cox[[8]]),
                                  HR = round(M.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(M.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(M.stage_sig_cox[[8]][,4],3),
                                  p.value=round(M.stage_sig_cox[[7]][,5],3))
M.stage_sig_cox_dat

#Stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat



#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     T.stage_sig_cox_dat,
                     N.stage_sig_cox_dat,
                     M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     # Grade_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "Gender",
                        "T.stage",
                        "N.stage",
                        "M.stage",
                        "Stage",
                        # "Grade",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('06_nomo/Fig6b.pdf',height = 4,width = 6,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='#009966',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 6,graph.pos = 2)
dev.off()

#
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.stage +  N.stage +M.stage+Stage+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c(
                         # "Gender",
                         "T.stage",
                         "N.stage",
                         "M.stage",
                         "Stage",
                         # "Grade",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('06_nomo/Fig6c.pdf',height = 4,width = 6,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='#009966',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 6,graph.pos = 2)
dev.off()

# pdf('06_nomo/Nomogram.pdf', width = 12, height = 10)
# nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
#                                 
#                                 T.stage=tcga_cox_datas$T.stage,
#                                 N.stage=tcga_cox_datas$N.stage,
#                                 
#                      os = tcga_cox_datas$OS.time,
#                      status = tcga_cox_datas$OS,
#                      mks = c(1,3,5)
# )
# dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))



####07IC50和gsea############
dir.create('07_IC50_and_GSEA')
##7.1 #######
load("07_IC50_and_GSEA/tcga_durg_ic50.RData")
tcga_durg_ic50_res[1:5,1:5]
colnames(tcga_durg_ic50_res)[1]='Cisplatin'

IC50.mat=data.frame(Riskscore=tcga.risktype.cli$Riskscore,tcga_durg_ic50_res[tcga.risktype.cli$Samples,])

library(ggcorrplot)
library(psych)
IC50_RS_cor <- corr.test(x =IC50.mat$Riskscore,
                         y = IC50.mat[,-1],
                         method = "spearman",adjust = "BH",ci = F)


IC50_RS_cor_res=data.frame(drugs=colnames( IC50.mat[,-1]))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
head(IC50_RS_cor_res)
table(IC50_RS_cor_res$p.adj<0.05,abs(IC50_RS_cor_res$cor)>0.2)
IC50_RS_cor_res=IC50_RS_cor_res[IC50_RS_cor_res$p.adj<0.05 & abs(IC50_RS_cor_res$cor)>0.2,]
IC50_RS_cor_res=IC50_RS_cor_res[order(IC50_RS_cor_res$cor),]
head(IC50_RS_cor_res)

library(rcartocolor)
pdf('IC50_RS_cor.pdf',height = 8,width = 6,onefile = F)
fig7a <- ggplot(data=IC50_RS_cor_res,aes(x=cor,y=reorder(drugs,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_colour_gradient(low ='#FFC100' ,high = "#C40C0C")+
  geom_segment(aes(yend=drugs,xend=0),size=.5) +
  labs(x='spearman Correlation',y='Drugs')+theme_bw()+
  theme(text = element_text(family = 'Times'))
dev.off()
#7.1.#####
##2.2 GSEA##############
tcga.geneList=getGeneFC(gene.exp=tcga.exp[,tcga.risktype.cli$Samples], group=tcga.risktype.cli$Risktype,ulab='High',dlab='Low')

h.all.gmt<-read.gmt("data/h.all.v2023.1.Hs.entrez.gmt")
set.seed(456)
tcga.geneList.gsea<-GSEA(tcga.geneList,TERM2GENE = h.all.gmt,seed=T)
gsea_res_symbol <- setReadable(tcga.geneList.gsea,"org.Hs.eg.db",keyType = "ENTREZID")


tcga.geneList.gsea.res=tcga.geneList.gsea@result
write.xlsx(tcga.geneList.gsea.res,'07_IC50_and_GSEA/TCGA.subtypeGSEA.res.xlsx',overwrite = T)
head(tcga.geneList.gsea.res)
table(tcga.geneList.gsea.res$NES>0)
tcga.geneList.gsea.res$group=ifelse(tcga.geneList.gsea.res$NES>0,"Activated","Suppressed")
topPathways = tcga.geneList.gsea.res %>% group_by(group) %>% slice_max(n =10, order_by = abs(NES))
# library(GseaVis)
#  
#  GseaVis::gseaNb(gsea_res_symbol, geneSetID  = 'HALLMARK_E2F_TARGETS')

source('data/my_gseaNb.R')
source('data/')
fig7b=my_mgseaNb(object = tcga.geneList.gsea,subPlot =2,
                 geneSetID = topPathways$ID[topPathways$NES>0],
                 termWidth = 35,kegg = F,
                 legend.position = c(0.8,0.8),
                 addPval = T,
                 pvalX = 0.05,pvalY = 0.05)
fig7c=my_mgseaNb(object = tcga.geneList.gsea,subPlot =3,
                 geneSetID = topPathways$ID[topPathways$NES<0],
                 termWidth = 35,kegg = F,
                 legend.position = c(0.8,0.8),
                 addPval = T,
                 pvalX = 0.05,pvalY = 0.05,
                 curveCol = c4a('brewer.paired',12))

pdf('07_IC50_and_GSEA/Fig7.pdf',height = 14,width = 16,onefile = F)
mg_merge_plot(fig7a,mg_merge_plot(fig7b,fig7c,labels = c('C','D'),nrow=2,heights = c(1,1.3)),
              ncol=2,widths =  c(1,1.5),labels = c('A',''))           

dev.off()
save.image("LUAD4.Rdata")
setwd('/home2/s421955/projects/retrovirus/data')
load('/home2/s421955/projects/retrovirus/data/all_mat_sato_cb.RData')
#load('/home2/s421955/projects/retrovirus/data/combinedsiggenes.RData')
load('~/temp/1107/cbhervgenes_p.RData')
hervst=ctrl_length+1
herved=dim(all_mat_cb)[1]
cancer_type='Clear cell'
samples_allfastq=samples_allfastq[samples_allfastq$Data_type!='DNA',]
#######(3.5) correlation between each genes
pdf('~/temp/0109/pairs_allcancer_sato.pdf')
herv=sig_all_p[[1]]$genes
if(length(herv)<2){next}
pairs(t(all_mat_cb[herv,]),pch=16,main=cancer_type)
dev.off()
##############cell infiltration level
cor_herv_sato=c()
cellinfilt=function(all_mat_cb,sig_herv){
  library("lme4")
  library("pbkrtest")
  load("~/projects/retrovirus/data/eTME_signatures_v4.RData")#change signatures to the newest eTME v4
  signatures=signatures[-21]#normal inputs involve no 'pDCs' class genes
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("GSVA")
  ssgsea_cb=function(exp_data,signatures){
    library("GSVA")
    ssgsea=sapply(names(signatures),function(cell)
      gsva(exp_data,list(signatures[[cell]]),method="ssgsea",rnaseq=T)[1,])
    return(ssgsea)
  }
  infilt_score=ssgsea_cb(all_mat_cb[1:ctrl_length,],signatures)
  #cor_herv=sapply(row.names(all_mat_cb)[hervst:herved],function(gene) cor(infilt_score,all_mat_cb[gene,],method = 'spearman'))
  cor_herv=cor(infilt_score,t(all_mat_cb[hervst:dim(all_mat_cb)[1],]),method = 'spearman')
  herv_combined_enrichment=gsva(data.matrix(all_mat_cb),list(sig_herv),method="ssgsea",rnaseq=T)[1,]
  cortest=cor(infilt_score,herv_combined_enrichment,method='spearman')
  names(cortest)=names(signatures)
  cortest=cortest[order(abs(cortest),decreasing = T)]
  row.names(cor_herv)=names(signatures)
  return(list(cor_herv,cortest))
}
if(grepl(cancer_type,'Oncocytoma')){sub_tumor=grepl(cancer_type,samples_allfastq$Histology)&grepl('^kidney tumor$',samples_allfastq$Source)
  }else{sub_tumor=grepl(paste('^',cancer_type,sep=''),samples_allfastq$Histology)&grepl('^kidney tumor$',samples_allfastq$Source)}
cor_herv_sato=cellinfilt(data.matrix(all_mat_cb[,sub_tumor]),sig_all_p[[1]]$genes)
plot_cor=function(corr,sig,cancer_type,celltype)
{
  corr=corr[celltype,]
  corr=corr[order(corr)]
  
  marked=list('cornflowerblue'=toupper(sig),'brown'=c('AF080231.1',"EU137846.2"))
  #marked=list('cornflowerblue'=c('AF080231.1',"EU137846.2",toupper(sig)))#This mark is for preslides, all interested hervs are in blue without words
  markedlabel=c('HERV-K','HERV-E')
  if(cancer_type=='FHD'){cancer_type='FHD/FHDL'}
  bp=barplot(corr,names="",col='gray88',border=NA,horiz=TRUE,main=cancer_type) 
  for (col in names(marked))
  {
    for (gene in marked[[col]]) 
    {
      if (!gene %in% names(corr)) {next}
      i=which(names(corr)==gene)
      y=bp[i,]
      segments(x0=0,x1=corr[i],y0=y,y1=y,col=col,lwd=1)
      if (col=="brown") {text(x=0,y=y,markedlabel[which(marked$brown==gene)],pos=2,cex=0.7)} # this is a bit ad hoc, but will work here
    }
  }
}
pdf('~/temp/0411/ssgsea_corr_sato_Th.pdf',width=4,height = 6)
cor2plot=cor_herv_sato[1]
plot_cor(cor2plot[[1]],sig_all_p[[1]]$genes,cancer_type,'Th cells')
dev.off()
############survival analysis
library(survival)
followup=read.csv(file='/project/bioinformatics/Xiao_lab/shared/genomics/Sato/survival.csv',row.names = 1,header = T,stringsAsFactors = F)
followup[followup[,2]=='alive',2]=0
followup[followup[,2]=='dead',2]=1
followup=data.matrix(followup[,c(2,3)])
patsurvid=row.names(followup)
patsurvid=sapply(patsurvid,function(name) gsub('ccRCC','Sato',name))
row.names(followup)=patsurvid
patlist_surv=intersect(row.names(followup),samples_allfastq$Patient_ID)
KCPsurvival=function(all_mat_cb,patlist_surv,samples_allfastq,cancer_type,followup,sig_all_p,hervst){
  if(cancer_type=='all'){
    sub_tumor=!grepl('kidney normal',samples_allfastq$Source)
  }else{
    if(length(samples_allfastq[samples_allfastq$Histology==cancer_type,1])==0){print('Invalid cancer type!');return(NA)
    }else{
      if(cancer_type=='Oncocytoma'){
        sub_tumor=grepl(cancer_type,samples_allfastq$Histology)&(!grepl('kidney normal',samples_allfastq$Source))
      }else{
        sub_tumor=grepl(paste('^',cancer_type,sep=''),samples_allfastq$Histology)&(!grepl('kidney normal',samples_allfastq$Source))
      }
    }
  }
  patlist_surv_tmp=intersect(unique(samples_allfastq[sub_tumor,"Patient_ID"]),patlist_surv)
  #patients have both abnormal samples and survival data in this cancer subtype
  if(length(patlist_surv_tmp)<30){print('Non-abundance')}
  all_mat_herv_sub=all_mat_cb[hervst:dim(all_mat_cb)[1],samples_allfastq[samples_allfastq$Patient_ID%in%patlist_surv_tmp,"ResultID"]]
  all_mat_herv_sub=aggregate(t(all_mat_herv_sub),by=list(samples_allfastq[samples_allfastq$Patient_ID%in%patlist_surv_tmp,"Patient_ID"]),mean)
  row.names(all_mat_herv_sub)=all_mat_herv_sub[,1]
  all_mat_herv_sub=data.matrix(all_mat_herv_sub[patlist_surv_tmp,-1])
  #Attention: always look at sequences when use %in%, two data may order differently.
  #remove patients with only tumor or nromal samples
  survtmp=Surv(followup[patlist_surv_tmp,2],followup[patlist_surv_tmp,1])
  if(all(table(followup[patlist_surv_tmp,2])==0)){print('No event.');next}
  #fit=survfit(survtmp~1,type='kaplan-meier')
  #plot(fit$time,fit$surv,type='s')
  hr=c()
  for(gene in sig_all_p$genes){
    fit=coxph(survtmp~all_mat_herv_sub[,gene]>median(all_mat_herv_sub[,gene]))
    pval=summary(fit)$coefficients[5]
    hrtmp=summary(fit)$conf.int[1]
    ci1=summary(fit)$conf.int[3]
    ci2=summary(fit)$conf.int[4]
    HRtmp=list(pval=pval,hrtmp=hrtmp,ci1=ci1,ci2=ci2)
    hr=cbind(hr,HRtmp)
  }
  colnames(hr)=sig_all_p$genes
  return(hr)
}
cox_sato=c()
all(row.names(patlist_surv)%in%samples_allfastq[samples_allfastq$Histology=='Clear cell',"Patient_ID"])
#all survival data belongs to ccRCC patients
cox_sato=KCPsurvival(all_mat_cb,patlist_surv,samples_allfastq,cancer_type,followup,sig_all_p[[1]],hervst)
save(cox_sato,file='~/temp/0109/survresult_sato.RData')
save(cor_herv_sato,file='~/temp/0109/corr_cell_sato.RData')

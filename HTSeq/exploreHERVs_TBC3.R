#Load the 2 hervselect functions from exploreHERVs_TBC2.R
setwd('/home2/s421955/projects/retrovirus/data')
#load('/home2/s421955/projects/retrovirus/data/all_mat_mergedherv.RData')
load('/home2/s421955/projects/retrovirus/data/all_mat_mergedherv_p.RData')#will use this piece for the whole process
load('/home2/s421955/projects/retrovirus/data/all_mat.RData')
hervst=ctrl_length+1
herved=dim(all_mat)[1]

#################compare paired_hervs between T and N
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
sig_all_p=c()
for(i in 1:length(cancer_type)){
  sig_all_p[[i]]=hervselect_paired(hervst,herved,cancer_type[i],all_mat,samples_allfastq,all_mat_cb_p,status_p,herv_group_p)
}
save(sig_all_p,file='~/temp/1107/cbhervgenes_p.RData')
load('~/temp/1107/noncbhervgenes_p.RData')
load('~/temp/1107/cbhervgenes_p.RData')
pdf(file='~/temp/0109/paired_herv_hox_cb.pdf',width=12,height=8.5)
par(mfrow=c(3,2))
for(i in 1:length(cancer_type)){
  type=sig_all_p[[i]]$cancer_type
  herv_sigtmp=sig_all_p[[i]]$genes
  patind=sig_all_p[[i]]$patients
  exp_tmpT=as.matrix(sig_all_p[[i]]$exp_hervT)
  exp_tmpN=as.matrix(sig_all_p[[i]]$exp_hervN)
  if(length(herv_sigtmp)==1){
    plot(x=c(1,1.5),y=c(exp_tmpN[1,1],exp_tmpT[1,1]),xlim=c(0.5,2),ylim=range(c(range(exp_tmpT),range(exp_tmpN))),main=type,xlab=NA,ylab='Expression Level',xaxt='n')
    axis(side=1,at=1.25,labels = NA)
    text(x=1.25,y=-1,srt=45,adj=1,xpd=TRUE,labels=herv_sigtmp)
    for(j in 1:length(patind)){
      points(x=1,y=exp_tmpN[j,1],pch=16,col='black')
      points(x=1.5,y=exp_tmpT[j,1],pch=16,col='red')
      segments(x0=1,x1=1.5,y0=exp_tmpN[j,1],y1=exp_tmpT[j,1],col='grey',lwd=2)
    }
  }else{
    if(type=='FHD'){type='FHD/FHDL'} #trick at here, need FHD in 'grepl' for include both FHD and FHDL, here change to mark title.
    plot(x=c(1,1.5),y=c(exp_tmpN[herv_sigtmp[1],1],exp_tmpT[herv_sigtmp[1],1]),xlim=c(0,length(herv_sigtmp)+1),ylim=c(0,12),main=type,xlab=NA,ylab='Expression Level',xaxt='n')
    axis(side=1,at=c(1:length(herv_sigtmp)+0.25),labels=NA)
    text(x=c(1:length(herv_sigtmp)+0.25),y=-1,srt=45,adj=1,xpd=TRUE,labels=herv_sigtmp)
    for(z in 1:length(herv_sigtmp)){
      for(j in 1:length(patind)){
        points(x=z,y=exp_tmpN[herv_sigtmp[z],j],pch=16,col='black')
        points(x=z+0.5,y=exp_tmpT[herv_sigtmp[z],j],pch=16,col='red')
        segments(x0=z,x1=z+0.5,y0=exp_tmpN[herv_sigtmp[z],j],y1=exp_tmpT[herv_sigtmp[z],j],col='grey',lwd=2)
      }
    }
  }
}
dev.off()
###############use paired hervs do unpaired bar plot########
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
sig_all=c()
for(i in 1:length(cancer_type)){
  sig_all[[i]]=hervselect_unpaired(hervst,herved,cancer_type[i],all_mat,samples_allfastq,all_mat_cb_p,status_p,herv_group_p,sig_all_p[[i]]$genes)
}#could use only paired test to select as much hervs as possible
pdf('~/temp/0109/herv_hox_cb.pdf',width=12,height=8.5)
par(mfrow=c(2,3))
for(i in 1:length(cancer_type)){
  type=sig_all[[i]]$cancer_type
  if(type=='FHD'){type='FHD/FHDL'}
  herv_sigtmp=sig_all[[i]]$genes
  mat2bar=data.frame(t(cbind(sig_all[[i]]$exp_hervT,sig_all[[i]]$exp_hervN)))
  mat2bar$type=factor(c(rep('Tumor',dim(sig_all[[i]]$exp_hervT)[2]),rep('Normal',dim(sig_all[[i]]$exp_hervN)[2])))
  mat2bar=melt(mat2bar)
  boxplot(value~type*variable,data=mat2bar,col=c('grey','red'),ylim=c(0,13),ylab='Expression Level',xaxt='n',main=type)
  #axis(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),las=3,labels = unique(mat2bar$variable))
  #mtext(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),text = unique(mat2bar$variable))
  axis(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),labels=NA)
  text(x=seq(2,2*length(unique(mat2bar$variable))+1,2),y=-1,srt=45,adj=1,xpd=TRUE,labels = unique(mat2bar$variable))
}
dev.off()
rm(mat2bar)

################heatmap showing ocmparasion among all sources
#sig_all don't contain herv_sig genes in FHD
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
sources=c('kidney normal','^kidney tumor$',' tumorgraft','metastasis$','^kidney thrombus$')
mat2heat=c()
for(i in 1:length(cancer_type)){
  mattemp=c()
  for(source in sources){
    matind=grepl(source,samples_allfastq$Source)&grepl(cancer_type[i],samples_allfastq$Histology)
    hervind=sig_all[[i]]$genes
    mattmp=all_mat_cb_p[hervind,matind]
    if(sum(matind)>1){
      mattmp=apply(mattmp,1,mean)
    }
    if(sum(matind)==0){
      mattmp=rep(0,length(hervind))
      print(paste('No samples:',cancer_type,source,sep=' '))
    }
    mattemp=cbind(mattemp,mattmp)
  }
  colnames(mattemp)=c('Normal','Tumor','Tumorgraft','Metastasis','Thrombus')
  if(cancer_type[i]=='FHD'){cancer_type[i]='FHD/FHDL'}
  row.names(mattemp)=paste(hervind,cancer_type[i],sep='|')
  mat2heat=rbind(mat2heat,mattemp)
}
#mat2heat=mat2heat[order(sapply(1:dim(mat2heat)[1],function(i) strsplit(row.names(mat2heat)[i],split='|',fixed = T)[[1]][1])),]
pdf('~/temp/0109/Allsourceheat_cb.pdf',height=8,width=7)
heatmap.2(data.matrix(mat2heat),margins=c(4,12),trace='none',col=c('black',brewer.pal(9,'RdBu')),breaks=c(0,0.5,1:9),notecex=0.5,Rowv=FALSE,Colv=FALSE,dendrogram ='none',srtCol=45,cexCol = 0.9,cexRow = 0.9)
dev.off()
#######(3.5) correlation between each genes
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
#sig_all don't contain herv_sig genes in FHD
#examine unpaired significant genes
sources=c('kidney normal','^kidney tumor$')
pdf('~/temp/0109/pairs_allcancer.pdf',width=20,height=20)
for(i in 1:length(cancer_type)){
  herv=sig_all_p[[i]]$genes
  if(length(herv)<2){next}
  exp_T=t(sig_all_p[[i]]$exp_hervT)
  exp_N=t(sig_all_p[[i]]$exp_hervN)
  mat2pairs=rbind(exp_N,exp_T)
  col=c(rep('black',dim(exp_N)[1]),rep('red',dim(exp_T)[1]))
  if(cancer_type[i]=='FHD'){cancer_type[i]='FHD/FHDL'}
  pairs(mat2pairs,col=col,pch=16,main=cancer_type[i])
} 
dev.off()
##############cell infiltration level
cor_herv=c()
cellinfilt=function(all_mat_cb_p,sig_herv){
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
  infilt_score=ssgsea_cb(all_mat_cb_p[1:ctrl_length,],signatures)
  #cor_herv=sapply(row.names(all_mat_cb_p)[hervst:herved],function(gene) cor(infilt_score,all_mat_cb_p[gene,],method = 'spearman'))
  cor_herv=cor(infilt_score,t(all_mat_cb_p[hervst:dim(all_mat_cb_p)[1],]),method = 'spearman')
  herv_combined_enrichment=gsva(data.matrix(all_mat_cb_p),list(sig_herv),method="ssgsea",rnaseq=T)[1,]
  cortest=cor(infilt_score,herv_combined_enrichment,method='spearman')
  names(cortest)=names(signatures)
  cortest=cortest[order(abs(cortest),decreasing = T)]
  row.names(cor_herv)=names(signatures)
  return(list(cor_herv,cortest))
}
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
for(i in 1:length(cancer_type)){
  if(grepl(cancer_type[i],'Oncocytoma')){sub_tumor=grepl(cancer_type[i],samples_allfastq$Histology)&grepl('^kidney tumor$',samples_allfastq$Source)
  }else{sub_tumor=grepl(paste('^',cancer_type[i],sep=''),samples_allfastq$Histology)&grepl('^kidney tumor$',samples_allfastq$Source)}
  cor_herv[[i]]=cellinfilt(data.matrix(all_mat_cb_p[,sub_tumor]),sig_all_p[[i]]$genes)
}
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
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
pdf('~/temp/0411/preslide_ssgsea_corr_herv_Th.pdf',width=4,height=6)
for(i in 1:length(cancer_type)){
  cor2plot=cor_herv[[i]]
  plot_cor(cor2plot[[1]],sig_all_p[[i]]$genes,cancer_type[i],'Th cells')
}
dev.off()
############survival analysis
library(survival)
IDcopath=read.csv(file='/home2/s421955/cleaned_data/sample.csv',row.names = 1,header = T,stringsAsFactors = F)
followup=read.csv(file='/home2/s421955/cleaned_data/newKCP-Followup.csv',row.names = 1,header = T,stringsAsFactors = F)
split_pat_id<-function(a,b)
{
  a=gsub("c[0-9]$","",a,perl=T)
  a=gsub("[abcdefg]","",a,perl=T)
  a=gsub('M[0-9]$','',a,perl=T)
  a=gsub("c1[0-9]$","",a,perl=T)
  a=toupper(a)
  b=toupper(b)
  c=strsplit(a,b)[[1]]
  paste(c,collapse=b)
}
ID2patient=sapply(1:dim(IDcopath)[1],function(i) split_pat_id(IDcopath$XPID[i],IDcopath$Type[i]))
patconvert=data.frame(patID=IDcopath$PatientID,trimmedID=ID2patient,stringsAsFactors = F)
patconvert=patconvert[!grepl('[_|TH|N]',patconvert$trimmedID,perl=T),]
patconvert=patconvert[!duplicated(patconvert$patID),]
patconvert$readID=NA
for(i in 1:dim(patconvert)[1]){
  if(!is.na(match(patconvert[i,2],samples_allfastq$Patient_ID))){patconvert$readID[i]=patconvert$trimmedID[i]}
  else{
    if(length(grep(patconvert[i,2],samples_allfastq$Root_tree))>0){
      patidtmp=samples_allfastq$Patient_ID[grep(patconvert[i,2],samples_allfastq$Root_tree)]
      if(length(unique(patidtmp))==1){patconvert$readID[i]=unique(patidtmp)}
      else{print(patconvert[i,])}
    }else{print(patconvert[i,])}
  }
}
patconvert$readID[patconvert$trimmedID=='219']='51'
patconvert$readID[patconvert$trimmedID=='189']='53'
patconvert$readID[patconvert$trimmedID=='310']='96'
patconvert=patconvert[!is.na(patconvert$readID),]
#patconvert matches between KCP IDs and Patient_ID in samples_allfastq
row.names(followup)=followup$PatientID

KCPsurvival=function(all_mat_cb_p,patconvert,samples_allfastq,cancer_type,followup,sig_all,hervst){
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
  patlist=intersect(unique(samples_allfastq[sub_tumor,"Patient_ID"]),patconvert$readID)
  #patients have both abnormal samples and survival data in this cancer subtype
  if(length(patlist)<30){print('Non-abundance')}
  sub_tumor_cleaned=sub_tumor&samples_allfastq$Patient_ID%in%patlist
  sub_survival_cleaned=patconvert$readID%in%patlist
  all_mat_herv_sub=all_mat_cb_p[hervst:dim(all_mat_cb_p)[1],sub_tumor_cleaned]
  all_mat_herv_sub=aggregate(t(all_mat_herv_sub),by=list(samples_allfastq$Patient_ID[sub_tumor_cleaned]),mean)
  row.names(all_mat_herv_sub)=all_mat_herv_sub[,1]
  all_mat_herv_sub=data.matrix(all_mat_herv_sub[patconvert[sub_survival_cleaned,"readID"],-1])
  #Attention: always look at sequences when use %in%, two data may order differently.
  #remove patients with only tumor or nromal samples
  survtmp=Surv(followup[patconvert[sub_survival_cleaned,"patID"],'DateLastFollowup'],followup[patconvert[sub_survival_cleaned,"patID"],"VitalStatus"])
  if(all(table(followup[patconvert[sub_survival_cleaned,"patID"],"VitalStatus"])==0)){print('No event.');next}
  ########temprary: draw ks surfit plots
  #group=c()
  #group=rep("high",dim(all_mat_herv_sub)[1])
  #group[all_mat_herv_sub[,gene]<=median(all_mat_herv_sub[,gene])]="low"
  #fit=survfit(survtmp~all_mat_herv_sub[,gene]>median(all_mat_herv_sub[,gene]),type='kaplan-meier')
  #fit=survfit(survtmp~group,type='kaplan-meier')
  #pdf('~/temp/0109/surv3_KCP_ccRCC.pdf',width=4,height = 4)
  #plot(fit,col=c('blue','orange'),conf.int=F,main=gene)
  #legend('topright',legend=c('high','low'),col=c('blue','orange'),lty=c(1,1),cex=0.8)
  #dev.off()
  hr=c()
  for(gene in sig_all$genes){
    fit=coxph(survtmp~all_mat_herv_sub[,gene]>median(all_mat_herv_sub[,gene]))
    pval=summary(fit)$coefficients[5]
    hrtmp=summary(fit)$conf.int[1]
    ci1=summary(fit)$conf.int[3]
    ci2=summary(fit)$conf.int[4]
    HRtmp=list(pval=pval,hrtmp=hrtmp,ci1=ci1,ci2=ci2)
    hr=cbind(hr,HRtmp)
  }
  colnames(hr)=sig_all$genes
  return(hr)
}
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
cox=c()
for(i in 1:length(cancer_type)){
  cox[[i]]=KCPsurvival(all_mat_cb_p,patconvert,samples_allfastq,cancer_type[i],followup,sig_all_p[[i]],hervst)
}
save(cox,file='~/temp/0109/survresult_KCP.RData')
save(cor_herv,file='~/temp/0109/corr_cell_KCP.RData')

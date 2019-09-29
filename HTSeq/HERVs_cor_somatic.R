#Feb 5th, 2018
#This script is to analyze the correlation between somatic mutations 
#include exonic and splicing mutations, exclude syn SNVs, genes focus on PBRM1, BAP1, mTOR1 involved in PD1 pathway
setwd('~/projects/retrovirus/data')
#source('~/PycharmProjects/HTSeq/exploreHERVs_TBC2.R')
#import hervselect_paired fucntions from source
load('~/temp/1107/noncbhervgenes_p.RData')
load('/home2/s421955/projects/retrovirus/data/all_mat_mergedherv_p.RData')#will use this piece for the whole process
load('/home2/s421955/projects/retrovirus/data/all_mat.RData')
#1. preprocess HERVs and exp
hervst=ctrl_length+1
herved=dim(all_mat)[1]
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
sig_all_paired_KCP=c()
for(i in 1:length(cancer_type)){
  sig_all_paired_KCP[[i]]=hervselect_paired(hervst,herved,cancer_type[i],all_mat,samples_allfastq,all_mat_cb_p,status_p,herv_group_p)
}
#2. preprocess somatic mutations for KCP
#/project/bioinformatics/Xiao_lab/shared/neoantigen/data/somatic/exome_seq
setwd('/project/bioinformatics/Xiao_lab/shared/neoantigen/data/somatic/exome_seq')
KCP_allmut=list()
KCP_mutfactmat=data.frame(matrix(NA,ncol=4,nrow=dim(samples_allfastq)[1]),stringsAsFactors = F)
colnames(KCP_mutfactmat)=c('Root_tree','PBRM1','BAP1','mTORc1')
KCP_mutfactmat$Root_tree=samples_allfastq[,"Root_tree"]
getMut=function(samples_allfastq,KCP_allmut,KCP_mutfactmat,column){
  #column='Root_tree' or 'PatientID', depends on cohort
  for(i in 1:dim(samples_allfastq)[1]){
  tmp=samples_allfastq[i,column]
  file_tmp=file.path(tmp,'somatic_mutations_hg38.txt')
  if(grepl('Sato',tmp)){file_tmp=file.path(tmp,'somatic_mutations_hg19.txt')}
  if(file.exists(file_tmp)){
    open_file=read.table(file_tmp,sep='\t',header=T,stringsAsFactors = F)
    ind1=grepl('^exonic',open_file$Func.refGene,perl=T)|grepl('splicing',open_file$Func.refGene)
    ind2=!grepl('^synonymous SNV',open_file$ExonicFunc.refGene)
    open_file=open_file[ind1&ind2,]
    KCP_allmut[[i]]=open_file
    names(KCP_allmut)[i]=tmp
    KCP_mutfactmat[i,1]=tmp
    KCP_mutfactmat[i,'PBRM1']=any(grepl('PBRM1',open_file$Gene.refGene,ignore.case = T))
    KCP_mutfactmat[i,'BAP1']=any(grepl('BAP1',open_file$Gene.refGene,ignore.case = T))
    KCP_mutfactmat[i,'mTORc1']=any(grepl('mTORc1',open_file$Gene.refGene,ignore.case = T))
    }else{print(paste('Unvalid:',tmp))}
  }
  return(list(mut=KCP_allmut,fact=KCP_mutfactmat))
}
result_KCP=getMut(samples_allfastq,KCP_allmut,KCP_mutfactmat,'Root_tree')
KCP_mutfactmat=na.omit(result_KCP$fact[,1:4])
KCP_allmut=result_KCP$mut[!sapply(result_KCP$mut,is.null)]
KCP_index=samples_allfastq
#manually remove one duplicated root_tree
KCP_index=KCP_index[-352,]
KCP_allmut=KCP_allmut[-124]
KCP_mutfactmat=KCP_mutfactmat[-124,]

#3. plot groups for KCP
#test:mutfactmat=KCP_mutfactmat;exp=all_mat_cb_p;sig=sig_all_paired_KCP;index=KCP_index;file='~/temp/testfig2';width=12;height=8
getplot=function(mutfactmat,exp,sig,index,file,width,height,column){
  pdf(file=file,width=width,height = height)
  mutfactmat$ind=sapply(mutfactmat[,column],function(ele) which(index[,column]==ele))
  for(i in 1:length(sig)){
    type=sig[[i]]$cancer_type
    matind=grepl(type,index[mutfactmat$ind,"Histology"])
    if(!sum(matind)>0){next}
    sig_genes=sig[[i]]$genes
    expind=index[mutfactmat$ind,"ResultID"][matind]
    expsub=exp[sig_genes,expind]
    mat2bar=cbind.data.frame(mutfactmat[matind,2:4],t(expsub))
    row.names(mat2bar)=mutfactmat[matind,1]
    mat2bar=melt(mat2bar)
    if(sum(mat2bar[,"PBRM1"])>0){
      boxplot(value~PBRM1*variable,data=mat2bar,col=c('grey85','indianred1'),ylab='Expression level',xaxt='n',main=paste(type,'PBRM1','mutation'))
      axis(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),labels=NA)
      text(x=seq(2,2*length(unique(mat2bar$variable))+1,2),y=-1,srt=45,adj=1,xpd=TRUE,labels = unique(mat2bar$variable))
    }
    if(sum(mat2bar[,"BAP1"])>0){
      boxplot(value~BAP1*variable,data=mat2bar,col=c('grey85','lightsalmon'),ylab='Expression level',xaxt='n',main=paste(type,'BAP1','mutation'))
      axis(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),labels=NA)
      text(x=seq(2,2*length(unique(mat2bar$variable))+1,2),y=-1,srt=45,adj=1,xpd=TRUE,labels = unique(mat2bar$variable))
    }
    if(sum(mat2bar[,"mTORc1"])>0){
      boxplot(value~mTPRc1*variable,data=mat2bar,col=c('grey85','palegreen'),ylab='Expression level',xaxt='n',main=paste(type,'mTORc1','mutation'))
      axis(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),labels=NA)
      text(x=seq(2,2*length(unique(mat2bar$variable))+1,2),y=-0.5,srt=45,adj=1,xpd=TRUE,labels = unique(mat2bar$variable))
    }
  }
  dev.off()
}
getplot(KCP_mutfactmat,all_mat_cb_p,sig_all_paired_KCP,KCP_index,file='~/temp/KCP_sommutgroup.pdf',width=5,height=5,column = 'Root_tree')

#4. preprocess for TCGA
load('/home2/s421955/projects/retrovirus/data/all_mat_TCGA_cb.RData')
TCGA_exp=all_mat_cb
TCGA_index=samples_allfastq[samples_allfastq$Data_type!='DNA'&samples_allfastq$Source=='kidney tumor',]
TCGA_index=TCGA_index[!duplicated(TCGA_index$Patient_ID),]
#TCGA cohort call somatic mutation with patient ID, use cancer samples to map HERV expression
hervst=ctrl_length+1
herved=dim(TCGA_exp)[1]
setwd('/project/bioinformatics/Xiao_lab/shared/neoantigen/data/somatic/exome_seq')
TCGA_allmut=list()
TCGA_mutfactmat=data.frame(matrix(NA,ncol=4,nrow=dim(TCGA_index)[1]),stringsAsFactors = F)
colnames(TCGA_mutfactmat)=c('Patient_ID','PBRM1','BAP1','mTORc1')
result_TCGA=getMut(TCGA_index,TCGA_allmut,TCGA_mutfactmat,'Patient_ID')
TCGA_mutfactmat=na.omit(result_TCGA$fact[,1:4])
TCGA_allmut=result_TCGA$mut[!sapply(result_TCGA$mut,is.null)]
getplot(TCGA_mutfactmat,TCGA_exp,sig_all_paired_KCP,TCGA_index,file='~/temp/TCGA_sommutgroup.pdf',width=5,height=5,column = 'Patient_ID')

#5. preprocess for Sato
load('/home2/s421955/projects/retrovirus/data/all_mat_sato_cb.RData')
Sato_exp=all_mat_cb
Sato_index=samples_allfastq[samples_allfastq$Data_type!='DNA'&samples_allfastq$Source=='kidney tumor',]
Sato_index=Sato_index[!duplicated(Sato_index$Patient_ID),]
hervst=ctrl_length+1
herved=dim(Sato_exp)[1]
setwd('/project/bioinformatics/Xiao_lab/shared/neoantigen/data/somatic/exome_seq')
Sato_allmut=list()
Sato_mutfactmat=data.frame(matrix(NA,ncol=4,nrow=dim(Sato_index)[1]),stringsAsFactors = F)
colnames(Sato_mutfactmat)=c('Patient_ID','PBRM1','BAP1','mTORc1')
result_Sato=getMut(Sato_index,Sato_allmut,Sato_mutfactmat,'Patient_ID')
Sato_mutfactmat=na.omit(result_Sato$fact[,1:4])
Sato_allmut=result_Sato$mut[!sapply(result_Sato$mut,is.null)]
getplot(Sato_mutfactmat,Sato_exp,sig_all_paired_KCP,Sato_index,file='~/temp/Sato_sommutgroup.pdf',width=5,height=5,column = 'Patient_ID')


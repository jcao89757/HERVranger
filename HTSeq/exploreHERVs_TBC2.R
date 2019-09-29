setwd('/home2/s421955/projects/retrovirus/data')
load('/home2/s421955/projects/retrovirus/data/all_mat.RData')#This contains data without combined but normalized
library(gplots)
library(RColorBrewer)
library(reshape2)
hervst=ctrl_length+1
herved=dim(all_mat)[1]
hervselect_paired=function(hervst,herved,cancer_type,exp_data,meta,exp_data_cb,status,herv_group_p){
  if(cancer_type=='all'){
    sub_tumor=grepl('^kidney tumor$',meta$Source)
    sub_normal=grepl('kidney normal',meta$Source)
  }else{
    if(length(meta[meta$Histology==cancer_type,1])==0){print('Invalid cancer type!');return(NA)
    }else{
      if(cancer_type=='Oncocytoma'){
        sub_tumor=grepl(cancer_type,meta$Histology)&grepl('^kidney tumor$',meta$Source)
        sub_normal=grepl(cancer_type,meta$Histology)&grepl('kidney normal',meta$Source)
      }else{
        sub_tumor=grepl(paste('^',cancer_type,sep=''),meta$Histology)&grepl('^kidney tumor$',meta$Source)
        sub_normal=grepl(paste('^',cancer_type,sep=''),meta$Histology)&grepl('kidney normal',meta$Source)
      }
    }
  }
  patlist=intersect(unique(meta[sub_tumor,"Patient_ID"]),unique(meta[sub_normal,'Patient_ID']))
  #patients have both tumor and normal samples in this cancer subtype
  sub_tumor_cleaned=sub_tumor&meta$Patient_ID%in%patlist
  #remove patients with only tumor or nromal samples
  sub_normal_cleaned=sub_normal&meta$Patient_ID%in%patlist
  herv_genes=row.names(exp_data)[hervst:herved]
  exp_mat_hervT=aggregate(t(exp_data[row.names(exp_data)%in%herv_genes,sub_tumor_cleaned]),by=list(meta$Patient_ID[sub_tumor_cleaned]),mean)
  rownames(exp_mat_hervT)=exp_mat_hervT[,1]
  exp_mat_hervT=t(as.matrix(exp_mat_hervT[,-1]))
  exp_mat_hervN=aggregate(t(exp_data[row.names(exp_data)%in%herv_genes,sub_normal_cleaned]),by=list(meta$Patient_ID[sub_normal_cleaned]),mean)
  rownames(exp_mat_hervN)=exp_mat_hervN[,1]
  exp_mat_hervN=t(as.matrix(exp_mat_hervN[,-1]))
  pvals=sapply(herv_genes,function(gene) t.test(exp_mat_hervT[gene,],exp_mat_hervN[gene,],alternative = 'greater',paired = TRUE)$p.val)
  folds=sapply(herv_genes,function(gene) {mean(exp_mat_hervT[gene,])-mean(exp_mat_hervN[gene,])})
  #original cutoff is folds>2 and pvals<0.001,11/07/2017
  sig_genes=herv_genes[folds>1.5&pvals<0.001]
  
  #dev.off()
  if(!is.na(exp_data_cb)){
    exp_mat_hervT=aggregate(t(exp_data_cb[row.names(exp_data_cb)%in%herv_genes,sub_tumor_cleaned]),by=list(meta$Patient_ID[sub_tumor_cleaned]),mean)
    rownames(exp_mat_hervT)=exp_mat_hervT[,1]
    exp_mat_hervT=t(as.matrix(exp_mat_hervT[,-1]))
    exp_mat_hervN=aggregate(t(exp_data_cb[row.names(exp_data_cb)%in%herv_genes,sub_normal_cleaned]),by=list(meta$Patient_ID[sub_normal_cleaned]),mean)
    rownames(exp_mat_hervN)=exp_mat_hervN[,1]
    exp_mat_hervN=t(as.matrix(exp_mat_hervN[,-1]))
    sig_genes=unique(sapply(status[[i]],function(ind) herv_group_p$herv_siggenes_p[herv_group_p$herv_sigstat_p==ind]))
  }
  sig=list(genes=sig_genes,cancer_type=cancer_type,patients=colnames(exp_mat_hervT),exp_hervT=exp_mat_hervT[sig_genes,],exp_hervN=exp_mat_hervN[sig_genes,])
  return(sig)
} 

hervselect_unpaired=function(hervst,herved,cancer_type,exp_data,meta,exp_data_cb,status,herv_group,presethervs=NA){
  if(cancer_type=='all'){
    sub_tumor=grepl('^kidney tumor$',meta$Source)
    sub_normal=grepl('kidney normal',meta$Source)
  }else{
    if(length(meta[meta$Histology==cancer_type,1])==0){print('Invalid cancer type!');return(NA)
    }else{
      if(cancer_type=='Oncocytoma'){
        sub_tumor=grepl(cancer_type,meta$Histology)&grepl('^kidney tumor$',meta$Source)
        sub_normal=grepl(cancer_type,meta$Histology)&grepl('kidney normal',meta$Source)
      }else{
        sub_tumor=grepl(paste('^',cancer_type,sep=''),meta$Histology)&grepl('^kidney tumor$',meta$Source)
        sub_normal=grepl(paste('^',cancer_type,sep=''),meta$Histology)&grepl('kidney normal',meta$Source)
      }
    }
  }
  herv_genes=row.names(exp_data)[hervst:herved]
  exp_mat_hervT=data.matrix(exp_data[herv_genes,sub_tumor])
  exp_mat_hervN=data.matrix(exp_data[herv_genes,sub_normal])
  pvals=sapply(herv_genes,function(gene) t.test(exp_mat_hervT[gene,],exp_mat_hervN[gene,],alternative = 'greater',paired=FALSE)$p.val)
  folds=sapply(herv_genes,function(gene) {mean(exp_mat_hervT[gene,])-mean(exp_mat_hervN[gene,])})
  #original cutoff is folds>2 and pvals<0.001,11/07/2017
  sig_genes=herv_genes[folds>1.5&pvals<0.001]
  if(!is.na(exp_data_cb)){
    
    exp_mat_hervT=aggregate(t(exp_data_cb[row.names(exp_data_cb)%in%herv_genes,sub_tumor]),by=list(meta$Patient_ID[sub_tumor]),mean)
    rownames(exp_mat_hervT)=exp_mat_hervT[,1]
    exp_mat_hervT=t(as.matrix(exp_mat_hervT[,-1]))
    exp_mat_hervN=aggregate(t(exp_data_cb[row.names(exp_data_cb)%in%herv_genes,sub_normal]),by=list(meta$Patient_ID[sub_normal]),mean)
    rownames(exp_mat_hervN)=exp_mat_hervN[,1]
    exp_mat_hervN=t(as.matrix(exp_mat_hervN[,-1]))
    sig_genes=unique(sapply(status[[i]],function(ind) herv_group$herv_siggenes[herv_group$herv_sigstat==ind]))
  }
  if(!is.na(presethervs)){sig_genes=presethervs;print('With presethervs');print(sig_genes)}
  if(length(sig_genes)>0){
    #sig=list(genes=sig_genes,cancer_type=cancer_type,exp_hervT=exp_data[sig_genes,sub_tumor],exp_hervN=exp_data[sig_genes,sub_normal])
    sig=list(genes=sig_genes,cancer_type=cancer_type,exp_hervT=exp_mat_hervT[sig_genes,],exp_hervN=exp_mat_hervN[sig_genes,])
  }else{sig=list(genes=sig_genes,cancer_type=cancer_type)}
  return(sig)
}
source=c('kidney normal','^kidney tumor$','^kidney tumorgraft$','metastasis','^kidney thrombus$')
sub_source=c('nromal','tumor','tumorgraft','metastasis','thrombus')
#Calculate sig_all and sig_all_p, with original herv annotation
cancer_type=c('Clear cell','Papillary','Chromophobe','Oncocytoma')
sig_all=c()
for(i in 1:length(cancer_type)){
  sig_all[[i]]=hervselect_unpaired(hervst,herved,cancer_type[i],all_mat,samples_allfastq,NA,NA)
}
save(sig_all,file='~/temp/1107/noncbhervgenes.RData')
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
sig_all_p=c()
for(i in 1:length(cancer_type)){
  sig_all_p[[i]]=hervselect_paired(hervst,herved,cancer_type[i],all_mat,samples_allfastq,NA)
}
#######correlation between each genes
# cancer_type=c('Clear cell','Papillary','Chromophobe','Oncocytoma')
# #sig_all don't contain herv_sig genes in FHD
# #examine unpaired significant genes
# sources=c('kidney normal','^kidney tumor$')
# pdf('~/temp/pairs_allcancer.pdf')
# for(i in 1:length(cancer_type)){
#   herv=sig_all[[i]]$genes
#   if(length(herv)<2){next}
#   exp_T=t(sig_all[[i]]$exp_hervT)
#   exp_N=t(sig_all[[i]]$exp_hervN)
#   mat2pairs=rbind(exp_N,exp_T)
#   col=c(rep('black',dim(exp_N)[1]),rep('red',dim(exp_T)[1]))
#   pairs(mat2pairs,col=col,pch=16,main=cancer_type[i])
# } 
# dev.off()
# status=list(c(1,1,1,1,2,3,3,4),c(5,6,7,7,8,9),c(10,11),c(12,12,10,13,13,13,13,12,12))
#observed from paired plots
#######correlation between each genes in paired search
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
sources=c('kidney normal','^kidney tumor$')
pdf('~/temp/1107/pairs_allcancer_p.pdf',height=30,width=30)
for(i in 1:length(cancer_type)){
  herv=sig_all_p[[i]]$genes
  if(length(herv)<2){next}
  exp_T=t(sig_all_p[[i]]$exp_hervT)
  exp_N=t(sig_all_p[[i]]$exp_hervN)
  mat2pairs=rbind(exp_N,exp_T)
  col=c(rep('black',dim(exp_N)[1]),rep('red',dim(exp_T)[1]))
  pairs(mat2pairs,col=col,pch=16,main=cancer_type[i])
} 
dev.off()
grouprf=read.table('~/temp/1107/hervcombine',stringsAsFactors = F)
#calc grouprf from seq2cluster.py in mymac/PycharmProjects/HTSeq/
status_p=c()
for(i in 1:length(sig_all_p)){
  status=sapply(sig_all_p[[i]]$genes,function(gene) grouprf[gene,1])
  status_p[[i]]=status
}
save(sig_all_p,status_p,file='~/temp/1107/noncbhervgenes_p.RData')
########combine hervs at here
load('/home2/s421955/projects/retrovirus/data/all_ctrl.RData')
load('/home2/s421955/projects/retrovirus/data/all_herv.RData')
annotategenes=function(exp_data){
  ann=read.table("~/data/genome/hg38/annotations.txt",stringsAsFactors = F)
  ann=ann[!duplicated(ann),]
  rownames(ann)=ann$V1
  exp_data=exp_data[rownames(exp_data) %in% ann$V1,]
  ann=ann[rownames(exp_data),]
  table(is.na(ann$V1))
  exp_data=aggregate(exp_data,by=list(ann$V2),sum) #sum up all isoforms of one gene
  rownames(exp_data)=exp_data[,1]
  exp_data=as.matrix(exp_data[,-1])
  return(exp_data)
}#Only ctrl genes can exist in ann!
ctrl_mat=annotategenes(ctrl_mat)
# herv_siggenes=unlist(sapply(sig_all,'[',1))
# herv_sigstat=unlist(status)
# herv_group=data.frame(herv_siggenes,herv_sigstat,stringsAsFactors = F)
# herv_group=herv_group[-19,]
# herv_sub2group=herv_mat[row.names(herv_mat)%in%herv_group[,1],]
# herv_subssub=herv_mat[!row.names(herv_mat)%in%herv_group[,1],]
# herv_sub2group=aggregate(herv_sub2group[herv_group[,1],],by=list(herv_group[,2]),sum)
# herv_sub2group=herv_sub2group[,-1]
# row.names(herv_sub2group)=herv_group[!duplicated(herv_group[,2]),1]
# herv_group=herv_group[!duplicated(herv_group[,2]),]
# herv_mat=rbind(herv_sub2group,herv_subssub)
# rm(herv_sub2group,herv_subssub)
# all_mat_cb=rbind(ctrl_mat,herv_mat)
# all_mat_cb[]=normalize.quantiles(log(data.matrix(all_mat_cb)+1))
# ctrl_length=dim(ctrl_mat)[1]
# save(all_mat_cb,ctrl_length,samples_allfastq,herv_group,file='/home2/s421955/projects/retrovirus/data/all_mat_mergedherv.RData')
#combine paired hervs
load('/home2/s421955/projects/retrovirus/data/all_herv.RData')
#Update at Jan 9th, 2018
#Use counts/length of seq (or exon seqs length normally) convert herv counts into herv counts/kb,to make different seqs compareble
cplength=read.table(file='cplength.txt',stringsAsFactors = F)
cpname=read.table(file='cpName.txt',stringsAsFactors = F)
row.names(cplength)=cpname[,1]
hervcp_p=cplength[row.names(herv_mat),1]
herv_mat[]=sapply(colnames(herv_mat),function(name) herv_mat[,name]/hervcp_p*1000)
#end of update
herv_siggenes_p=unlist(sapply(sig_all_p,'[',1))
herv_sigstat_p=unlist(status_p)
herv_group_p=data.frame(herv_siggenes_p,herv_sigstat_p,stringsAsFactors = F)
herv_group_p=herv_group_p[!duplicated(herv_group_p[,1]),]
herv_sub2group_p=herv_mat[row.names(herv_mat)%in%herv_group_p[,1],]
herv_subssub_p=herv_mat[!row.names(herv_mat)%in%herv_group_p[,1],]
herv_sub2group_p=aggregate(herv_sub2group_p[herv_group_p[,1],],by=list(herv_group_p[,2]),sum)
herv_sub2group_p=herv_sub2group_p[,-1]
row.names(herv_sub2group_p)=herv_group_p[!duplicated(herv_group_p[,2]),1]
herv_group_p=herv_group_p[!duplicated(herv_group_p[,2]),]
herv_mat=rbind(herv_sub2group_p,herv_subssub_p)
all_mat_cb_p=rbind(ctrl_mat,herv_mat)
all_mat_cb_p[]=normalize.quantiles(log(data.matrix(all_mat_cb_p)+1))
ctrl_length=dim(ctrl_mat)[1]
rm(ctrl_mat,herv_mat)
rm(herv_sub2group,herv_subsub,herv_sub2group_p,herv_subssub_p)

save(all_mat_cb_p,ctrl_length,samples_allfastq,herv_group_p,file='/home2/s421955/projects/retrovirus/data/all_mat_mergedherv_p.RData')
############produce plots in 
############Unpaired hervs between T and N
############Don't have to clear before TBC3.R
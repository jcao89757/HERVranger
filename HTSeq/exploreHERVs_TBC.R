#This script is to further analyze HERV exp_data
#--outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 20
#######(0) HERV CNV or percentage? 
# (1) global analysis, significant HERVs,unpaired T and N, in ccRCC. 
# (2) Apply them / find sth new in T and N paired (within same group of patients)
# (3) expand to other histology and source. 
# (3.5) All significant HERVs together and investigate correlations of expression within them across histologies. 
# (4) campare differnet histologies (5) correlate with survival 
# (6) correlate with infiltration (7) validate in Sato and TCGA 
setwd('/home2/s421955/projects/retrovirus/data')
load('/home2/s421955/projects/retrovirus/data/raw_all_mat.RData')#all_mat without log and normalization
load('/home2/s421955/projects/retrovirus/data/all_mat_mergedherv_test.RData')#all_mat merged high correlated herv
load('/home2/s421955/projects/retrovirus/data/all_mat.RData')#all_mat normalized with original hervs
library(gplots)
library(RColorBrewer)
hervst=ctrl_length+1
herved=dim(all_mat)[1]
hervselect_paired=function(hervst,herved,cancer_type,exp_data,meta){
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
  exp_mat_hervT=aggregate(t(exp_data[herv_genes,sub_tumor_cleaned]),by=list(meta$Patient_ID[sub_tumor_cleaned]),mean)
  rownames(exp_mat_hervT)=exp_mat_hervT[,1]
  exp_mat_hervT=t(as.matrix(exp_mat_hervT[,-1]))
  exp_mat_hervN=aggregate(t(exp_data[herv_genes,sub_normal_cleaned]),by=list(meta$Patient_ID[sub_normal_cleaned]),mean)
  rownames(exp_mat_hervN)=exp_mat_hervN[,1]
  exp_mat_hervN=t(as.matrix(exp_mat_hervN[,-1]))
  pvals=sapply(herv_genes,function(gene) t.test(exp_mat_hervT[gene,],exp_mat_hervN[gene,],alternative = 'greater',paired = TRUE)$p.val)
  folds=sapply(herv_genes,function(gene) {mean(exp_mat_hervT[gene,])-mean(exp_mat_hervN[gene,])})
  #change at here: original cutoff is folds>2 and pvals<0.001,11/07/2017
  sig_genes=herv_genes[folds>1.5&pvals<0.001]
  #dev.off()
  sig=list(genes=sig_genes,cancer_type=cancer_type,patients=colnames(exp_mat_hervT),exp_hervT=exp_mat_hervT[sig_genes,],exp_hervN=exp_mat_hervN[sig_genes,])
  return(sig)
} 

hervselect_unpaired=function(hervst,herved,cancer_type,exp_data,meta){
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
  #change at here: original cutoff is folds>2 and pvals<0.001,11/07/2017
  sig_genes=herv_genes[folds>1.5&pvals<0.001]
  if(length(sig_genes)>0){
    sig=list(genes=sig_genes,cancer_type=cancer_type,exp_hervT=exp_data[sig_genes,sub_tumor],exp_hervN=exp_data[sig_genes,sub_normal])
  }else{sig=list(genes=sig_genes,cancer_type=cancer_type)}
  return(sig)
}

source=c('kidney normal','^kidney tumor$','^kidney tumorgraft$','metastasis','^kidney thrombus$')
sub_source=c('nromal','tumor','tumorgraft','metastasis','thrombus')


#####(1) unpaired all type cancers * hervs in all types of cancers, box plot, only N and T
cancer_type=c('Clear cell','Papillary','Chromophobe','Oncocytoma')
#temprarely remove FHD because we didn;'t find sig hervs in this histology
#still no FHD sig_hervs after merge hervs
sig_all=c()
for(i in 1:length(cancer_type)){
  sig_all[[i]]=hervselect_unpaired(hervst,herved,cancer_type[i],all_mat,samples_allfastq)
}
pdf('~/temp/1107/unpaired_herv_hox_noncb.pdf')
#unpaired_herv_hox.pdf represents results without merging hervs in probably the same group
par(mfrow=c(2,2))
for(i in 1:length(cancer_type)){
  type=sig_all[[i]]$cancer_type
  herv_sigtmp=sig_all[[i]]$genes
  mat2bar=data.frame(t(cbind(sig_all[[i]]$exp_hervT,sig_all[[i]]$exp_hervN)))
  mat2bar$type=factor(c(rep('Tumor',dim(sig_all[[i]]$exp_hervT)[2]),rep('Normal',dim(sig_all[[i]]$exp_hervN)[2])))
  mat2bar=melt(mat2bar)
  boxplot(value~type*variable,data=mat2bar,col=c('grey','red'),ylab='Expression Level',xaxt='n',main=type)
  #axis(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),las=3,labels = unique(mat2bar$variable))
  #mtext(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),text = unique(mat2bar$variable))
  axis(side=1,at=seq(1.5,2*length(unique(mat2bar$variable))+0.5,2),labels=NA)
  text(x=seq(2,2*length(unique(mat2bar$variable))+1,2),y=-1,srt=45,adj=1,xpd=TRUE,labels = unique(mat2bar$variable))
}
dev.off()
rm(mat2bar)

######(2)Paired all hervs* cancer types, patients 
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')
sig_all_p=c()
for(i in 1:length(cancer_type)){
  sig_all_p[[i]]=hervselect_paired(hervst,herved,cancer_type[i],all_mat,samples_allfastq)
}
pdf('~/temp/paired_herv_hox_new.pdf')
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
    plot(x=c(1,1.5),y=c(exp_tmpN[herv_sigtmp[1],1],exp_tmpT[herv_sigtmp[1],1]),xlim=c(0,length(herv_sigtmp)+1),ylim=range(c(range(exp_tmpT),range(exp_tmpN))),main=type,xlab=NA,ylab='Expression Level',xaxt='n')
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
######(3.5)
cancer_type=c('Clear cell','Papillary','FHD','Chromophobe','Oncocytoma')

all_sig=c();all_norm=c()
for(i in 1:length(cancer_type)){all_sig=union(all_sig,sig_all[[i]]$genes);all_norm=union(all_norm,colnames(sig_all[[i]]$exp_hervN))}
mat2pairs=matrix(NA,ncol=6,nrow = length(all_sig))
for(i in 1:length(cancer_type)){
  patlist=colnames(sig_all[[i]]$exp_hervT)
  mat2pairs[,i]=apply(all_mat[all_sig,patlist],1,mean)
}
colnames(mat2pairs)=c(cancer_type,'Normal')
mat2pairs[,length(cancer_type)+1]=apply(all_mat[all_sig,all_norm],1,mean)
row.names(mat2pairs)=all_sig
#manually assign exp_tumors to FHD
sub_tumor=grepl('FHD',samples_allfastq$Histology)&grepl('^kidney tumor$',meta$Source)
mat2pairs[,3]=apply(all_mat[all_sig,sub_tumor],1,mean)
pairs(mat2pairs)#How to lable them?
mat2pairspar=matrix(NA,ncol=5,nrow = length(all_sig))
for(i in 1:length(cancer_type)){mat2pairspar[,i]=mat2pairs[,i]-mat2pairs[,6]}
colnames(mat2pairspar)=colnames(mat2pairs)[1:5]
row.names(mat2pairspar)=row.names(mat2pairs)
MASS::parcoord(mat2pairspar,col='grey',var.label = TRUE,lwd=2)
###########(3)
cancer_type=c('Clear cell','Papillary','Chromophobe','Oncocytoma')
#sig_all don't contain herv_sig genes in FHD
sources=c('kidney normal','^kidney tumor$',' tumorgraft','metastasis$','^kidney thrombus$')
mat2heat=c()
for(i in 1:length(cancer_type)){
  mattemp=c()
  for(source in sources){
    matind=grepl(source,samples_allfastq$Source)&grepl(cancer_type[i],samples_allfastq$Histology)
    hervind=sig_all[[i]]$genes
    mattmp=all_mat[hervind,matind]
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
  row.names(mattemp)=paste(hervind,cancer_type[i],sep='|')
  mat2heat=rbind(mat2heat,mattemp)
}
#mat2heat=mat2heat[order(sapply(1:dim(mat2heat)[1],function(i) strsplit(row.names(mat2heat)[i],split='|',fixed = T)[[1]][1])),]
pdf('~/temp/Allsourceheat_new.pdf',height=7,width=7)
heatmap.2(data.matrix(mat2heat),margins=c(4,12),trace='none',col=c('black',brewer.pal(9,'RdBu')),breaks=c(0,0.5,1:9),notecex=0.5,Rowv=FALSE,Colv=FALSE,dendrogram ='none',srtCol=45,cexCol = 0.9,cexRow = 0.9)
dev.off()
#######(3.5) correlation between each genes
cancer_type=c('Clear cell','Papillary','Chromophobe','Oncocytoma')
#sig_all don't contain herv_sig genes in FHD
#examine unpaired significant genes
sources=c('kidney normal','^kidney tumor$')
pdf('~/temp/pairs_allcancer.pdf')
for(i in 1:length(cancer_type)){
  herv=sig_all[[i]]$genes
  if(length(herv)<2){next}
  exp_T=t(sig_all[[i]]$exp_hervT)
  exp_N=t(sig_all[[i]]$exp_hervN)
  mat2pairs=rbind(exp_N,exp_T)
  col=c(rep('black',dim(exp_N)[1]),rep('red',dim(exp_T)[1]))
  pairs(mat2pairs,col=col,pch=16,main=cancer_type[i])
} 
dev.off()
#combine hervs based on pairs plot
sig_all[[1]]$status=c(1,1,1,1,2,3,3,4)
sig_all[[2]]$status=c(1,2,3,3,4,5)
sig_all[[3]]$status=c(1,2)
sig_all[[4]]$status=c(1,1,2,3,3,3,3,1,1)

#correlation test
#Currentlly won't use this as a combine reference, but rather pairs plot
# cors=cor(t(all_mat[hervst:herved,]),method='spearman')
# groups=rep(0,dim(cors)[1])
# label=0
# for(i in 1:dim(cors)[1]){
#   if(groups[i]==0){ 
#     label=label+1
#     groups[i]=label
#   }
#   for(j in i:dim(cors)[1])
#   {
#     if(i==j||groups[i]==groups[j]){next}
#     if(cors[i,j]>0.94){
#       if(groups[j]==0) {groups[j]=groups[i]} else {groups[groups==groups[j]]=groups[i]}
#     }
#   }
# }
# grouplist=aggregate(rownames(all_mat[hervst:herved,]),by=list(groups),function(x) paste(x,collapse=","))
# row.names(grouplist)=grouplist[,1]
# grouplist=grouplist[,-1]
# save(grouplist,file='~/projects/retrovirus/data/grouphervs.RData')

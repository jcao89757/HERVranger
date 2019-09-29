# you learned really well, these codes look quite good now, save for some minor improvements to be made
load('/home2/s421955/projects/retrovirus/data/all_mat.RData')
Nlist_fastq=samples_allfastq[samples_allfastq$Source=='kidney normal',]$ResultID # include normal of other tissue
Tlist_fastq=samples_allfastq[samples_allfastq$Source=='kidney tumor',]$ResultID 
TGlist_fastq=samples_allfastq[grepl('tumorgraft',samples_allfastq$Source),]$ResultID
Mlist_fastq=samples_allfastq[grepl('metastasis',samples_allfastq$Source),]$ResultID
Mlist_fastq=setdiff(Mlist_fastq,TGlist_fastq) # Metastasis tumorgraft should be included in tumorgraft part
###########Nlist_fastq contains different types of RNA sample file names########
all_mat_N=all_mat[,samples_allfastq$Source=='kidney normal']
all_mat_T=all_mat[,samples_allfastq$Source=='kidney tumor']
all_mat_TG=all_mat[,samples_allfastq$ResultID%in%TGlist_fastq]
all_mat_M=all_mat[,samples_allfastq$ResultID%in%Mlist_fastq]
###########all_mat_N contains different types of RNA counts quantile-normalized results############# 
###########all_mat_N contains different types of RNA counts quantile-normalized results############# 
library(gplots)
library(RColorBrewer)
################ttest betweein overall Normal and Tumor samples, hervst:herved indicates hervs##########
hervst=ctrl_length+1
herved=dim(all_mat)[1]
#################log scale for ttest#################
ttest_NT=sapply(hervst:herved,function(i) t.test(all_mat_N[i,],all_mat_T[i,],alternative='less',paired=F)$p.value)
sig_herv_ttest=as.vector(row.names(all_mat)[hervst:herved])[which(ttest_NT<0.01)]
###############T/N ratio of all patients average counts larger than 
all_mat_N_ave=apply(all_mat_N,1,mean)
all_mat_T_ave=apply(all_mat_T,1,mean)
ratioT2N=all_mat_T_ave-all_mat_N_ave
ratioT2N_herv=ratioT2N[hervst:herved]
herv_genes=row.names(all_mat)[hervst:herved]
herv_bycount=herv_genes[ratioT2N_herv>1]
herv_selected=intersect(herv_bycount,sig_herv_ttest)
all_mat_N_selected=data.matrix(all_mat_N[row.names(all_mat)%in%herv_selected,])
all_mat_T_selected=data.matrix(all_mat_T[row.names(all_mat)%in%herv_selected,])
all_mat_selected=t(cbind(all_mat_N_selected,all_mat_T_selected))
rm(all_mat_N,all_mat_T,all_mat_TG,all_mat_M)
#heatmap for all samples of selected hervs and each herv's density plot in all samples
#Not that useful.

#label=c(rep('Normal',dim(all_mat_N_selected)[2]),rep('Tumor',dim(all_mat_T_selected)[2]))
#label_color=rgb_palette(2)[as.factor(label)]
#pdf("~/temp/sig_herv.pdf")
#rgb_palette = colorRampPalette(c("green","black","red"), space = "rgb")
#heatmap.2(all_mat_selected,col=colorRampPalette(c("red", "yellow", "green"))(n = 299),key=TRUE,keysize=1.5,
  #scale='none',trace="none",dendrogram='column',Rowv=F,RowSideColors=label_color)
#dev.off()
#pdf('~/temp/dist.pdf')
#par(mfrow=c(2,13))
#for(i in 1:length(herv_selected)){
  #plot(density(all_mat_N_selected[i,]),type='l',lwd=1,col='black',main=herv_selected[i])
  #lines(density(all_mat_T_selected[i,]),type='l',lwd=1,col='red')
  #legend('topleft',legend=c('Normal','Tumor'),col=c('Black','Red'),lty=c(1,1),cex=0.5)
#}
#dev.off()
########################
hervselect=function(hervst,herved,cancer_type,exp_data,meta){
  #test:cancer_type='Clear cell',exp_data=all_mat,meta=samples_allfastq,
  sub_tumor=meta$Histology==cancer_type&(grepl('kidney tumor',meta$Source))
  #index of meta within both certain cancer_type and tumor
  sub_normal=meta$Histology==cancer_type&(grepl('kidney normal',meta$Source))
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
  rownames(exp_mat_hervN)=exp_mat_hervN[1]
  exp_mat_hervN=t(as.matrix(exp_mat_hervN[,-1]))
  pvals=sapply(herv_genes,function(gene) t.test(exp_mat_hervT[gene,],exp_mat_hervN[gene,],alternative = 'greater')$p.val)
  folds=sapply(herv_genes,function(gene) {mean(exp_mat_hervT[gene,])-mean(exp_mat_hervN[gene,])})
  sig_genes=herv_genes[folds>2&pvals<0.001]
  mat2heatmap=cbind(data.matrix(exp_data[sig_genes,sub_cleaned_meta$Result_T]),data.matrix(exp_data[sig_genes,sub_cleaned_meta$Result_N]))
  patcolor=c(rep('red',length(sub_cleaned_meta$Result_T)),rep('black',length(sub_cleaned_meta$Result_N)))
  #pdf(paste('~/temp/',gsub(' ','',cancer_type),'_heat.pdf',sep=''))
  heatmap.2(mat2heatmap,trace='none',col=brewer.pal(11,'RdBu'),ColSideColors=patcolor,Rowv=TRUE,Colv=FALSE,dendrogram ='row',key=TRUE,keysize=1.5,main=cancer_type)
  #dev.off()
  sig=list(genes=sig_genes,patients=row.names(sub_cleaned_meta))
  return(sig)
} 
cancer_type=list('Clear cell','Papillary','FHD','Chromophobe','Oncocyt')
sig=list()
for(i in 1:length(cancer_type)){
  sig[i]=hervselect(hervst,herved,cancer_type[i],all_mat,samples_allfastq)
}
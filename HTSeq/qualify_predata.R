#03202018
#This script is for the preliminary data of qualifying exam
setwd('/home2/s421955/projects/retrovirus/data')
load('/home2/s421955/projects/retrovirus/data/all_mat_mergedherv_p.RData')
library("lme4")
library("pbkrtest")
library('reshape')
#0. select 40 patients, each with 1 tumor sample
ind_tmp=(!duplicated(samples_allfastq$Patient_ID))&(grepl('tumor',samples_allfastq$Source))&(samples_allfastq$Histology=='Clear cell')
samples_allfastq=samples_allfastq[ind_tmp,]
#1. calculate their neoantigen and somatic mutation length
samples_allfastq$neo_num=NA
samples_allfastq$somatic_num=NA
path_neo='/project/bioinformatics/Xiao_lab/shared/neoantigen/data/neoantigen'
path_som='/project/bioinformatics/Xiao_lab/shared/neoantigen/data/somatic/exome_seq/'
name_neo='neoantigen_final.txt'
name_som='somatic_mutations_hg38.txt'
for(i in 1:dim(samples_allfastq)[1]){
  print(i)
  tmp=samples_allfastq$Root_tree[i]
  tmp_som_path=file.path(path_som,tmp,name_som)
  tmp_neo_path=file.path(path_neo,tmp,name_neo)
  if(file.exists(tmp_som_path)&file.exists(tmp_neo_path)){
    som_file=read.csv(tmp_som_path,stringsAsFactors = F,header = F)
    neo_file=read.csv(tmp_neo_path,stringsAsFactors = F,header = T)
    samples_allfastq$somatic_num[i]=dim(som_file)[1]-1
    samples_allfastq$neo_num[i]=dim(neo_file)[1]
  }else{print(tmp)}
}
samples_allfastq=na.omit(samples_allfastq)
corr=cor.test(samples_allfastq$neo_num,samples_allfastq$somatic_num,method='pearson')$estimate
corr=round(corr,2)
pdf(file='~/temp/somaticnut_neoantigen_cor.pdf',height=4,width=5)
plot(log(samples_allfastq$somatic_num,base=10)~c(1:dim(samples_allfastq)[1]),type='p',pch=16,col='lightsalmon',xaxt='n',main=NA,xlab=NA,ylab='log(n)',ylim=c(0,4))
points(log(samples_allfastq$neo_num,base=10)~c(1:dim(samples_allfastq)[1]),pch=16,col='navy')
axis(side=1,at=c(1:dim(samples_allfastq)[1]),labels = NA)
text(x=c(1:dim(samples_allfastq)[1]),y=-0.6,srt=45,adj=1,xpd=TRUE,labels=samples_allfastq$Patient_ID,cex=0.6)
text(x=5,y=0.5,labels=paste('Pearson correlation=',corr),cex=0.7)
legend('topright',legend=c('somatic mutation counts','neoantigen counts'),col=c('lightsalmon','navy'),pch=c(16,16),bty='n',cex=0.6)
dev.off()
#2. barplot of immune score, showing neoantigens are processed and immune active
load('/home2/s421955/data/all_data.rdata')
all_id=unlist(sapply(all_data,'[',7))
all_neo_imscore=matrix(NA,ncol=dim(samples_allfastq)[1],nrow=729)
colnames(all_neo_imscore)=samples_allfastq$Patient_ID
for(i in 1:dim(samples_allfastq)[1]){
  idtmp=samples_allfastq$Root_tree[i]
  if(length(which(all_id==idtmp))>0){
    ind_id=which(all_id==idtmp)
    imscore_tmp=all_data[[ind_id]]$neoantigen$immune
    all_neo_imscore[1:length(imscore_tmp),i]=imscore_tmp
  }
}
all_neo_imscore=all_neo_imscore[,-c(1,9)]
labels=colnames(all_neo_imscore)
all_neo_imscore=melt(all_neo_imscore[])
pdf(file='~/temp/neoantigen_immuscore.pdf',height=3,width=4)
boxplot(value~X2,data=all_neo_imscore,col='lightsalmon',ylim=c(-0.4,0.4),ylab='immunogenicity scores',xaxt='n')
text(x=c(1:length(labels)),y=-0.45,srt=45,adj=1,xpd=TRUE,labels=labels,cex=0.6)
abline(h=0,lty=2,col='navy')
text(x=2.2,y=-0.37,labels = 'non-epitopes',col='navy',cex=0.6)
abline(h=0.1,lty=2,col='navy')
text(x=1.5,y=0.35,labels = 'epitopes',col='navy',cex=0.6)
dev.off()
#3. ssgsea, infiltration of B/T cells
all_mat=data.matrix(all_mat_cb_p[,samples_allfastq$ResultID])
load("~/projects/retrovirus/data/eTME_signatures_v4.RData")
signatures_forall=c()
for(i in 1:length(signatures)){signatures_forall[[i]]=signatures[[i]][signatures[[i]]%in%row.names(all_mat)]}
names(signatures_forall)=names(signatures)
signatures_forall=signatures_forall[-21]
ssgsea_cb=function(exp_data,signatures){
  library("GSVA")
  ssgsea=sapply(names(signatures),function(cell)
    gsva(exp_data,list(signatures[[cell]]),method="ssgsea",rnaseq=T)[1,])
  return(ssgsea)
}
infilt_score=ssgsea_cb(all_mat,signatures_forall)
cor_neo_infilt_score=sapply(1:dim(infilt_score)[2],function(i) cor.test(samples_allfastq$neo_num,infilt_score[,i],method='spearman')$estimate)
cor_som_infilt_score=sapply(1:dim(infilt_score)[2],function(i) cor.test(samples_allfastq$somatic_num,infilt_score[,i],method='spearman')$estimate)
mat2dot=data.frame(neo=cor_neo_infilt_score,som=cor_som_infilt_score,stringsAsFactors = F)
row.names(mat2dot)=colnames(infilt_score)
pdf(file='~/temp/neoantigen_infilt_cor.pdf',height=5,width=4)
dotchart(x=mat2dot$neo,pch=16,labels=row.names(mat2dot),col='navy',main='Neoantigens',xlim=c(-0.5,0.2),xlab='spearman rho')
abline(v=0,lty=2,col='navy')
dev.off()
pdf(file='~/temp/somaticmutation_infilt_cor.pdf',height=5,width=4)
dotchart(x=mat2dot$som,pch=16,labels=row.names(mat2dot),col='lightsalmon',main='Somatic mutations',xlim=c(-0.5,0.2),xlab='spearman rho')
abline(v=0,lty=2,col='lightsalmon')
dev.off()
#4. signatures from Winslow
signatures_win=read.csv('signatures_winslow.csv',header=T,stringsAsFactors = F,quote = '')
signatures_win_agg=as.list(signatures_win)
signatures_win_agg$`fibroblasts`=toupper(c(signatures_win[,"X1"],signatures_win[,"X2"]))
signatures_win_agg$`endothelium`=toupper(c(signatures_win[,"X4"],signatures_win[,"X5"]))
signatures_win_agg$`immune`=toupper(as.character(unlist(c(signatures_win[,c(3,6:16)]))))
for(i in 1:length(signatures_win_agg)){signatures_win_agg[[i]]=signatures_win_agg[[i]][signatures_win_agg[[i]]%in%row.names(all_mat)]}
signatures_win_agg=signatures_win_agg[-5]
infilt_score_win=ssgsea_cb(all_mat,signatures_win_agg)
cor_neo_infilt_score=sapply(1:dim(infilt_score_win)[2],function(i) cor.test(samples_allfastq$neo_num,infilt_score_win[,i],method='spearman')$estimate)
cor_som_infilt_score=sapply(1:dim(infilt_score_win)[2],function(i) cor.test(samples_allfastq$somatic_num,infilt_score_win[,i],method='spearman')$estimate)
mat2dot=data.frame(neo=cor_neo_infilt_score,som=cor_som_infilt_score,stringsAsFactors = F)
row.names(mat2dot)=colnames(infilt_score_win)
pdf(file='~/temp/neoantigen_infilt_cor.pdf',height=5,width=4)
dotchart(x=mat2dot$neo,pch=16,labels=row.names(mat2dot),col='navy',main='Neoantigens',xlim=c(-0.6,0.2),xlab='spearman rho')
abline(v=0,lty=2,col='navy')
dev.off()
pdf(file='~/temp/somaticmutation_infilt_cor.pdf',height=5,width=4)
dotchart(x=mat2dot$som,pch=16,labels=row.names(mat2dot),col='lightsalmon',main='Somatic mutations',xlim=c(-0.6,0.2),xlab='spearman rho')
abline(v=0,lty=2,col='lightsalmon')
dev.off()
#5. survival plots
setwd('/home2/s421955/projects/retrovirus/data')
#load('/home2/s421955/projects/retrovirus/data/all_mat_mergedherv.RData')
load('/home2/s421955/projects/retrovirus/data/all_mat_mergedherv_p.RData')#will use this piece for the whole process
hervst=ctrl_length+1
herved=dim(all_mat_cb_p)[1]
load('~/temp/1107/cbhervgenes_p.RData')
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
cancer_type='Clear cell'
sub_tumor=grepl(paste('^',cancer_type,sep=''),samples_allfastq$Histology)&(!grepl('kidney normal',samples_allfastq$Source))
patlist=intersect(unique(samples_allfastq[sub_tumor,"Patient_ID"]),patconvert$readID)
sub_tumor_cleaned=sub_tumor&samples_allfastq$Patient_ID%in%patlist
sub_survival_cleaned=patconvert$readID%in%patlist
all_mat_herv_sub=all_mat_cb_p[hervst:dim(all_mat_cb_p)[1],sub_tumor_cleaned]
all_mat_herv_sub=aggregate(t(all_mat_herv_sub),by=list(samples_allfastq$Patient_ID[sub_tumor_cleaned]),mean)
row.names(all_mat_herv_sub)=all_mat_herv_sub[,1]
all_mat_herv_sub=data.matrix(all_mat_herv_sub[patconvert[sub_survival_cleaned,"readID"],-1])
row.names(followup)=followup$PatientID
survtmp=Surv(followup[patconvert[sub_survival_cleaned,"patID"],'DateLastFollowup'],followup[patconvert[sub_survival_cleaned,"patID"],"VitalStatus"])
if(all(table(followup[patconvert[sub_survival_cleaned,"patID"],"VitalStatus"])==0)){print('No event.');next}
group=c()
group=rep("high",dim(all_mat_herv_sub)[1])
gene="EU137846.2"
group[all_mat_herv_sub[,gene]<=median(all_mat_herv_sub[,gene])]="low"
fit=survfit(survtmp~all_mat_herv_sub[,gene]>median(all_mat_herv_sub[,gene]),type='kaplan-meier')
pdf('~/temp/0109/surv_HERVE_KCP_ccRCC.pdf',width=4,height = 4)
plot(fit,col=c('blue','orange'),conf.int=F,main='Survival plot for HERV-E',xlab='Days from last followup',ylab='Survival rate',cex=0.8)
legend('topright',legend=c('low','high'),col=c('blue','orange'),lty=c(1,1),cex=0.8)
dev.off()


#######add ResultID to samples_allfastq########
#######Realignment_qnorm.R --> exploreHERVs.R#############
setwd('/home2/s421955/projects/retrovirus/data')
#source("https://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
library('preprocessCore')

samples_allfastq=read.csv('/home2/s421955/cleaned_data/samples_allfastq.csv',stringsAsFactors = F,row.names = 1)
#samples_allfastq=read.csv('/home2/s421955/cleaned_data/samples_test.csv',stringsAsFactors = F,row.names = 1)#test for sato and TCGA
#FHD and FHDL, two resembel kinds of kidney cancer
RNAseqlist=sapply(1:dim(samples_allfastq)[1],function(i) paste(samples_allfastq$Root_tree[i],samples_allfastq$Seq_ID[i],sep='_'))
samples_allfastq$ResultID=RNAseqlist
samples_allfastq=samples_allfastq[!duplicated(samples_allfastq$ResultID),]

########save samples_allfastq as csv########

########Combine existed count results#######
data_path='/home2/s421955/projects/retrovirus/data'
ctrl='unmated_realignmentAligned.out.sam.ctrlfreqs.csv'
hervs='unmated_realignmentAligned.out.sam.featurefreqs.csv'
ctrl_mat=data.frame(matrix(ncol=1,nrow = 0))
herv_mat=data.frame(matrix(ncol=1,nrow = 0))
for(i in 1:dim(samples_allfastq)[1]){
  RNAseq=samples_allfastq$ResultID[i]
  ctrl_path=paste(data_path,'/',RNAseq,'/',ctrl,sep='')
  herv_path=paste(data_path,'/',RNAseq,'/',hervs,sep='')
  STARpath=paste(data_path,RNAseq,'/','unmated_realignment_STARtmp',sep='')
   if(file.exists(ctrl_path)&file.exists(herv_path)){
     ctrl_file=read.csv(ctrl_path,sep=',',header=T)
     herv_file=read.csv(herv_path,sep=',',header=T)
     colnames(ctrl_file)=c(NA,RNAseqlist[i])
     colnames(herv_file)=c(NA,RNAseqlist[i])
     # Ze, please make sure to talk to me tomorrow about this part. I want to double check
     # I guess you deleted some rows in the featureCounts output file. 
     # That is not necessary and could lead to some serious errors if not handled correctly
     # but it is good to see that you can work with more "advanced" commands in R now
     ctrl_mat=merge(ctrl_mat,ctrl_file,by.x=1,by.y=1,all.x=T,all.y=T,sort=T)
     herv_mat=merge(herv_mat,herv_file,by.x=1,by.y=1,all.x=T,all.y=T,sort=T)
   }else{print(RNAseq)} # good habit, always check for error and report loudly 
 }
##########combined results with rownames as geneIDs, colnames as ResultsIDs. 
row.names(ctrl_mat)=as.vector(ctrl_mat[,1])
ctrl_mat=data.matrix(ctrl_mat[,-1])
ctrl_mat[is.na(ctrl_mat)]=0 # sure?
row.names(herv_mat)=as.vector(herv_mat[,1])
herv_mat=herv_mat[,-1]
herv_mat[is.na(herv_mat)]=0
#########finished combination, annotate row names########
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
ctrl_mat_ann=annotategenes(ctrl_mat)
all_mat=rbind(ctrl_mat_ann,herv_mat)
save(all_mat,file='~/projects/retrovirus/data/raw_all_mat.RData')
all_mat[]=normalize.quantiles(log(data.matrix(all_mat)+1))
ctrl_length=dim(ctrl_mat_ann)[1]
save(ctrl_mat,file='/home2/s421955/projects/retrovirus/data/all_ctrl.RData')
save(herv_mat,file='/home2/s421955/projects/retrovirus/data/all_herv.RData')
save(samples_allfastq,all_mat,ctrl_length,file='/home2/s421955/projects/retrovirus/data/all_mat.RData')

#10/25/17 changes made at here.
#aggregate hervs by correlation (grouplist.RData)
load('/home2/s421955/projects/retrovirus/data/all_ctrl.RData')
load('/home2/s421955/projects/retrovirus/data/all_herv.RData')
load('~/projects/retrovirus/data/grouphervs.RData')
ctrl_mat=annotategenes(ctrl_mat)
grouplist=strsplit(grouplist,',')
status=sapply(1:length(grouplist),function(i) rep(i,length(grouplist[[i]])))
grouplist=unlist(grouplist)
status=unlist(status)
herv_mat=herv_mat[grouplist,]
herv_mat=aggregate(herv_mat,by=list(status),sum)
row.names(herv_mat)=grouplist[!duplicated(status)]
herv_mat=herv_mat[,-1]
#####aggregate all hervs in one group
all_mat=rbind(ctrl_mat,herv_mat)
all_mat[]=normalize.quantiles(log(data.matrix(all_mat)+1))
ctrl_length=dim(ctrl_mat)[1]
save(all_mat,ctrl_length,samples_allfastq,file='/home2/s421955/projects/retrovirus/data/all_mat_mergedherv.RData')
rm(ctrl_mat,ctrl_mat_ann,herv_mat)
#############after 10/05/2017 all_mat.RData become normalized, annotated RData
###############data preprocess end at here###########################
###############data analysis part moves to ~/PycharmProjects/HTSeq/exploreHERVs.R#### 



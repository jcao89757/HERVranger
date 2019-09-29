#######add ResultID to samples_allfastq########
#######Realignment_qnorm.R --> exploreHERVs.R#############
setwd('/home2/s421955/projects/retrovirus/data')
load('~/temp/1107/noncbhervgenes_p.RData')
#source("https://bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
library('preprocessCore')
samples_allfastq=read.csv('/home2/s421955/cleaned_data/samples_Sato_v2.csv',stringsAsFactors = F)#test for sato and TCGA
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
save(all_mat,file='~/projects/retrovirus/data/raw_all_mat_sato.RData')
all_mat[]=normalize.quantiles(log(data.matrix(all_mat)+1))
ctrl_length=dim(ctrl_mat_ann)[1]
save(samples_allfastq,all_mat,ctrl_length,file='/home2/s421955/projects/retrovirus/data/all_mat_sato_noncb.RData')

#10/25/17 changes made at here.
#aggregate hervs by correlation (grouplist.RData)

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
rm(herv_sub2group_p,herv_subssub_p)
all_mat_cb=rbind(ctrl_mat_ann,herv_mat)
ctrl_length=dim(ctrl_mat_ann)[1]
all_mat_cb[]=normalize.quantiles(log(data.matrix(all_mat_cb)+1))
save(all_mat_cb,ctrl_length,samples_allfastq,file='/home2/s421955/projects/retrovirus/data/all_mat_sato_cb.RData')
rm(ctrl_mat,ctrl_mat_ann,herv_mat)
#############after 10/05/2017 all_mat.RData become normalized, annotated RData
###############data preprocess end at here###########################
###############data analysis part moves to ~/PycharmProjects/HTSeq/exploreHERVs[0-9].R#### 



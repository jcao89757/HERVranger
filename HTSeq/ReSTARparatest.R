#This script is to test a proper pair of parameters for 2nd STAR alignment
#--outFilterMismatchNoverLmax 0.05/0.1/0.3 --outFilterMultimapNmax 5/20/50
#53 patient's normals and ccRCC tumor samples, total 106, are to be used in this test.
#make sure EU137846.2/EU137847.1 have the largest fold change in tumor/normal
library('preprocessCore')
setwd('/home2/s421955/projects/retrovirus/data')
meta_test=read.csv('/home2/s421955/cleaned_data/samples2test1007.csv',header=T,row.names=1,stringsAsFactors = F)

preprocess=function(ctrl,herv,meta){
  ctrl_mat=data.frame(matrix(ncol=1,nrow = 0))
  herv_mat=data.frame(matrix(ncol=1,nrow=0))
  for(i in 1:dim(meta)[1]){
    ctrl_file=paste(meta$ResultID[i],'_test/',ctrl,sep='')
    herv_file=paste(meta$ResultID[i],'_test/',herv,sep='')
    if(file.exists(ctrl_file)&file.exists(herv_file)){
      ctrl_exp=read.csv(ctrl_file,sep=',',header=T)
      herv_exp=read.csv(herv_file,sep=',',header=T)
      colnames(ctrl_exp)=c(NA,meta[i,'ResultID'])
      colnames(herv_exp)=c(NA,meta[i,'ResultID'])
      ctrl_mat=merge(ctrl_mat,ctrl_exp,by.x=1,by.y=1,all.x=T,all.y=T,sort=T)
      herv_mat=merge(herv_mat,herv_exp,by.x=1,by.y=1,all.x=T,all.y=T,sort=T)
    }else{print(paste('Unfound csv files:',meta[i,"ResultID"]))}
  }
  row.names(ctrl_mat)=ctrl_mat[,1]
  ctrl_mat=data.matrix(ctrl_mat[,-1])
  ctrl_mat[is.na(ctrl_mat)]=0 
  row.names(herv_mat)=herv_mat[,1]
  herv_mat=data.matrix(herv_mat[,-1])
  herv_mat[is.na(herv_mat)]=0
  all_mat=rbind(ctrl_mat,herv_mat)
  all_mat[]=normalize.quantiles(log(data.matrix(all_mat)+1))
  ctrl_length=dim(ctrl_mat)[1]
  print(paste('Combine finished:',ctrl,sep=''))
  resulttmp=list(exp_data=all_mat,ctrl_genes_num=ctrl_length,herv_file=herv)
  return(resulttmp)
}
countfolds=function(all_mat,meta,herv){
  if(is.element('EU137846.2',row.names(all_mat))){
    pval846=t.test(all_mat['EU137846.2',meta$Source=='kidney tumor'],all_mat['EU137846.2',meta$Source=='kidney normal'],alternative = 'greater',paired = T)$p.val
    fold846=mean(all_mat['EU137846.2',meta$Source=='kidney tumor'])-mean(all_mat['EU137846.2',meta$Source=='kidney normal'])
  }else{pval846=NA;fold846=NA;print(paste('No EU37846.2 counted',herv))}
  
  if(is.element('EU137847.1',row.names(all_mat))){
    pval847=t.test(all_mat['EU137847.1',meta$Source=='kidney tumor'],all_mat['EU137847.1',meta$Source=='kidney normal'],alternative = 'greater',paired = T)$p.val
    fold847=mean(all_mat['EU137847.1',meta$Source=='kidney tumor'])-mean(all_mat['EU137847.1',meta$Source=='kidney normal'])
  }else{pval847=NA;fold847=NA;print(paste('No EU37846.2 counted',herv))}
  resulttmp2=list(pval846=pval846,pval847=pval847,fold846=fold846,fold847=fold847)
  return(resulttmp2)
}

exp_all=c()
n=0
for(i in c(5,20,50)){
  for (j in c(0.05,0.1,0.3)){
    n=n+1
    ctrl_tmp=paste('unmated_realigntest_',i,'_',j,'Aligned.out.sam.ctrlfreqs.csv',sep='')
    herv_tmp=paste('unmated_realigntest_',i,'_',j,'Aligned.out.sam.featurefreqs.csv',sep='')
    #print(ctrl_tmp)
    exp_tmp=preprocess(ctrl_tmp,herv_tmp,meta_test)
    exp_tmp2=countfolds(exp_tmp$exp_data,meta_test,herv_tmp)
    exp_alltmp=list(exp_tmp,exp_tmp2)
    exp_all[[n]]=exp_alltmp
  }
}
save(exp_all,file='~/temp/tmp.RData')
load('~/temp/tmp.RData')

for(i in 1:9){#This part can display pvals and folds of all parameter pairs.
  print(exp_all[[i]][[1]][3])
  print(exp_all[[i]][[2]])
}
#exp_all contains 9 lists based on 9 parameter condations.
#exp_all[[i]] i=1:9,
#exp_all[[i]][[1]][1]=preprocessed exp_data,
#exp_all[[i]][[1]][2]=length of ctrl genes,
#exp_all[[i]][[1]][3]=one file name to look at the parameters' value
#exp_all[[i]][[2]]=pvals and folds for EU***846/847

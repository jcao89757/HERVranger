library(forestplot)
setwd('/home2/s421955/projects/retrovirus/data')
load('~/temp/0109/survresult_TCGA.RData')
load('~/temp/0109/survresult_sato.RData')
load('~/temp/0109/survresult_KCP.RData')
cox=cox[[1]]#will only compare within ccRCC patients among 3 cohorts
cox=matrix(unlist(t(cox)),ncol=4,byrow=FALSE,dimnames=list(colnames(cox),row.names(cox)))
pdf('~/temp/0109/forestsurvivalKCP.pdf',height=6,width=5)
forestplot(labeltext=cbind(c('HERVs','\n',row.names(cox)),c('P-values','\n',round(cox[,1],digits=3))),graph.pos=2,mean=c(NA,NA,cox[,2]),lower=c(NA,NA,cox[,3]),upper=c(NA,NA,cox[,4]),
           title='KCP ',xlab='<-Beneficial- ------Harmful------------>           ',
           txt_gp = fpTxtGp(label = gpar(cex=0.8),ticks=gpar(cex=0.8),xlab=gpar(cex=0.8),title=gpar(cex=1)),
           col=fpColors(box='black',lines='grey',zero='red'),
           zero=1,cex=0.9,lineheight='auto',boxsize=0.2,lwd.ci=2,
           ci.vertices=TRUE,ci.vertices.height = 0.4)
dev.off()

cox_sato=matrix(unlist(t(cox_sato)),ncol=4,byrow=FALSE,dimnames=list(colnames(cox_sato),row.names(cox_sato)))
pdf('~/temp/0109/forestsurvivalSato.pdf',height=6,width=5)
forestplot(labeltext=cbind(c('HERVs','\n',row.names(cox_sato)),c('P-values','\n',round(cox_sato[,1],digits=3))),graph.pos=2,mean=c(NA,NA,cox_sato[,2]),lower=c(NA,NA,cox_sato[,3]),upper=c(NA,NA,cox_sato[,4]),
           title='Sato',xlab='<-Beneficial- ------Harmful------------>           ',
           txt_gp = fpTxtGp(label = gpar(cex=0.8),ticks=gpar(cex=0.8),xlab=gpar(cex=0.8),title=gpar(cex=1)),
           col=fpColors(box='black',lines='grey',zero='red'),
           zero=1,cex=0.9,lineheight='auto',boxsize=0.2,lwd.ci=2,
           ci.vertices=TRUE,ci.vertices.height = 0.4)
dev.off()

cox_TCGA=matrix(unlist(t(cox_TCGA)),ncol=4,byrow=FALSE,dimnames=list(colnames(cox_TCGA),row.names(cox_TCGA)))
pdf('~/temp/0109/forestsurvivalTCGA.pdf',height=6,width=5)
forestplot(labeltext=cbind(c('HERVs','\n',row.names(cox_TCGA)),c('P-values','\n',round(cox_TCGA[,1],digits=3))),graph.pos=2,mean=c(NA,NA,cox_TCGA[,2]),lower=c(NA,NA,cox_TCGA[,3]),upper=c(NA,NA,cox_TCGA[,4]),
           title='TCGA',xlab='<-Beneficial- ------Harmful------------>           ',
           txt_gp = fpTxtGp(label = gpar(cex=0.8),ticks=gpar(cex=0.8),xlab=gpar(cex=0.8),title=gpar(cex=1)),
           col=fpColors(box='black',lines='grey',zero='red'),
           zero=1,cex=0.9,lineheight='auto',boxsize=0.2,lwd.ci=2,
           ci.vertices=TRUE,ci.vertices.height = 0.4)
dev.off()

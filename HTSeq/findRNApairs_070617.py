 #! usr/bin/env python
# This scratch is to find paired end RNA seqs from /project/BICF/shared/Kidney/rna/RAW
# The RNA IDs are read from sample.csv
#Scratch revised on 05/11/2017 to fit STAR as alignment tools

import pandas as pd
import numpy as np
import os

samples=pd.read_csv('/home2/s421955/cleaned_data/sample.csv',sep=',',header=0,names=None,dtype=str)
RNAseqlist=pd.read_csv('/home2/s421955/projects/retrovirus/data/index.txt',header=None,names=None,squeeze=True)#should contain R1 and R2
sample_pairs_R1=pd.Series()
sample_pairs_R2=pd.Series()
for i in range(0,len(samples['RNAID'])):
    RNAIDtmp=samples['RNAID'][i]
    if isinstance(RNAIDtmp,float):
        print('error')
        continue
    #read one row of RNAID each time
    RNAtmp=RNAseqlist[RNAseqlist.str.contains(RNAIDtmp)]
    #select R1+R2
    RNAtmpR1=RNAtmp[RNAtmp.str.contains('R1')].sort_values()
    RNAtmpR2 = RNAtmp[RNAtmp.str.contains('R2')].sort_values()
    if len(RNAtmpR1)!=len(RNAtmpR2):
        continue
    #some RNAIDs may contain several R1R2 pairs
    sample_pairs_R1=sample_pairs_R1.append(RNAtmpR1)
    sample_pairs_R2=sample_pairs_R2.append(RNAtmpR2)
sample_pairs_R1.index=range(0,len(sample_pairs_R1))
sample_pairs_R2.index=range(0,len(sample_pairs_R2))
iternum=0
#re-index, otherwise cannot access right item through index number
for i in range(0,len(sample_pairs_R1)):
    pathR1='/archive/BICF/shared/Kidney/rna/RAW/'+sample_pairs_R1[i]
    pathR2 = '/archive/BICF/shared/Kidney/rna/RAW/' + sample_pairs_R2[i]
    output='/home2/s421955/projects/retrovirus/data/'+sample_pairs_R1[i]+'.SAM'
    indexname='/home2/s421955/data/genome/hg38mm10/STAR'
    cmd='STAR --runMode alignReads --runThreadN 40 --genomeDir '+indexname+' --readFilesIn '+pathR1+' '+pathR2+\
        ' --readFilesCommand zcat --outFilterMultimapNmax 200 --outSAMunmapped Within --outSAMorder PairedKeepInputOrder --outReadsUnmapped None --outFileNamePrefix '+\
        '/home2/s421955/projects/retrovirus/data/'+sample_pairs_R1[i]+'/'+sample_pairs_R1[i]+'\n'
    filename=sample_pairs_R1[i]+'_070617.sh'
    if iternum==0:
        newfile=open(filename,'w')
        newfile.write('#!/bin/bash \n #SBATCH --job-name=serial \n #SBATCH --partition=super \n #SBATCH --nodes=1 \n #SBATCH --time=10-00:00:00 \n #SBATCH --output=./sbatch_output_%j\n')
    if iternum<40:
        newfile.write(cmd)
        iternum=iternum+1
        continue
    newfile.close()
    iternum=0
    #os.system('sbatch --partition=super '+filename)
    #os.system('srun '+filename)#not in need


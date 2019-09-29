#! usr/bin/env python
#This file is to (1) write shell cmd for realignment
#(2) continue if part of the jobs failed
#(3) if insert bam files, convert to fastq first, without gzip
import os
import sys
import glob
import pandas as pd
import numpy as np
index_file=sys.argv[1]
#index_file='/home2/s421955/cleaned_data/samples_genentech_nucleus.csv'
index_done=sys.argv[2]
#index_done='/home2/s421955/projects/retrovirus/data/index_done.txt'
#null done: '/home2/s421955/cleaned_data/fake.csv' for jobs start from the beginning
data_path=sys.argv[3]
#data_path='/home2/s421955/projects/retrovirus/data'
#output data path
index=pd.read_table(index_file,header=0,sep=',')
index_do=pd.read_table(index_done,header=None,sep=',')
j=0
for i in range(0,len(index)):
    line=index.iloc[i,:]
    if line['Data_type'].find('RNA')<0:
        continue
    if pd.isnull(line['Seq_ID']):
        continue
    indtmp=index_do[0].str.contains(line['Seq_ID'])
    filename=index_do[0][indtmp]
    newfilename='_'.join([line['Root_tree'],line['Seq_ID']])
    if not filename.empty:
        filenamevalue=filename.values
        if os.path.exists(os.path.join(data_path,filenamevalue[0])):
            os.rename(os.path.join(data_path,filenamevalue[0]),os.path.join(data_path,newfilename))
    else:
        nameclust=glob.glob(os.path.join(line['Path'],('*'+line['Seq_ID']+'*')))
        if len(nameclust)>2:
            print ('Multiple target RNA seq files: '+line['Seq_ID'])
            break
        elif len(nameclust)==0:
            print('Not found: '+SeqIDtmp)
            break
        else:
            pass
        for names in nameclust:
            if names.find('R1')>0 or names.find('_1_1')>0:
                name_R1=names
            if names.find('R2')>0 or names.find('_1_2')>0:
                name_R2=names
        cmd='python Realignment_main.py /home2/s421955/data/genome/hg38mm10/STAR ' \
            '/home2/s421955/data/genome/hg38mm10/STAR_HERV_092717 /home2/s421955/data/genome/hg38mm10/hg38mm10.gtf /home2/s421955/projects/retrovirus/data '+name_R1+' '+name_R2+' '+line['Path']+' '+newfilename+'_test'+' test'
        #print(cmd)
        if j==0:
            shfilename='Th_'+newfilename+'_cmd.sh'
            newfile=open(shfilename,'w')
            newfile.write(
              '#!/bin/bash \n #SBATCH --job-name=serial \n #SBATCH --partition=super \n #SBATCH --nodes=1 \n #SBATCH --time=10-00:00:00 \n #SBATCH --output=./sbatch_output_%j\n')
        if j<3:
            newfile.write(cmd)
            newfile.write('\n')
            j=j+1
        if j==3 or i==len(index)-1:
            newfile.close()
            j=0





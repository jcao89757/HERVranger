#! usr/bin/env python
# This is file is to build shell scripts for Realignment cmds.
import re
index=open('/Users/zhangze/Downloads/index.txt','r')
lines=index.readlines()
iternum=0
i=0
for name in lines:
    print(name)
    name_R1=re.split('[\t]',name)[0]
    name_R2=name_R1.replace('R1','R2')
    filename=name_R1+'.sh'
    if iternum%4==0:
        newfile = open(filename, 'w')
        newfile.write(
            '#!/bin/bash \n #SBATCH --job-name=serial \n #SBATCH --partition=super \n #SBATCH --nodes=1 \n #SBATCH --time=10-00:00:00 \n #SBATCH --output=./sbatch_output_%j\n')
    cmd='python Realignment_main.py /home2/s421955/data/genome/hg38mm10/STAR ' \
        '/home2/s421955/data/genome/hg38mm10/STAR_HERV_061517 /home2/s421955/data/genome/hg38mm10/hg38mm10.gtf ' \
        '/home2/s421955/projects/retrovirus/data '+name_R1+' '+name_R2+' /archive/BICF/shared/Kidney/rna/RAW'
    newfile.write(cmd)
    print(cmd)
    i=i+1
    newfile.write('\n')
    iternum = iternum + 1
    if iternum%4==0:
        newfile.close()
    #iternum = 0
index.close()

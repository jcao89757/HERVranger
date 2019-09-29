#! usr/bin/env python
# This is file is to build shell scripts for Realignment cmds.
index=open('/home2/s421955/projects/retrovirus/data/index.txt','r')
lines=index.readlines()
iternum=0
for name in lines:
    name = name.split('\n')[0]
    filename=name+'.sh'
    if iternum==0:
        newfile = open(filename, 'w')
        newfile.write(
            '#!/bin/bash \n #SBATCH --job-name=serial \n #SBATCH --partition=super \n #SBATCH --nodes=1 \n #SBATCH --time=10-00:00:00 \n #SBATCH --output=./sbatch_output_%j\n')

    if iternum<40:
        AlignRef='/home2/s421955/data/genome/hg38mm10/STAR_HERV_061517'
        Data_path='/home2/s421955/projects/retrovirus/data'
        Fasta_path = '/archive/BICF/shared/Kidney/rna/RAW'
        cmd_python=' '.join(['python Realignment_061617.py',AlignRef,Data_path,name,Fasta_path,'\n'])
        newfile.write(cmd_python)
        iternum = iternum + 1
        continue
    newfile.close()
    iternum = 0
index.close()

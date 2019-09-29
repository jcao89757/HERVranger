#! usr/bin/env python
# Ze_052217 This scratch's purpose includes (1)read sam files from STAR results with index.txt. (2) Mapped SAM files to
# gtf annotations. (3) Generate multimapped sequences from SAM files. (4) Generate unmapped short sequences from
# gzUnmapped files. (5) Wrap all unmaped files into .fq, with the names from index.txt.

import os
import numpy as np
index=open('/home2/s421955/projects/retrovirus/data/index.txt','r')
lines=index.readlines()
iternum=0
#use annotation downloaded from NCBI RefSeq
#GTF has removed HERVs
an_nonHERV='/home2/s421955/data/genome/hg38mm10/hg38mm10.gtf'

for name in lines:
    name=name.split('\n')[0]
    SAMfile='/home2/s421955/projects/retrovirus/data/'+name+'/'+name+'Aligned.out.sam'
    featureout= '/home2/s421955/projects/retrovirus/data/'+name +'/'+name+'Features.txt'
    feature_details= SAMfile+'.featureCounts'

    # -p: paired end default GTF annotation -M count Multimapped reads in BAM
    # -s 1/2 strand-specific read counting -T threads 1-32
     # count features
    cmd = 'featureCounts -T 32 -p -f --donotsort -R -M -a ' + an_nonHERV + ' -o ' + featureout + ' ' + SAMfile + '\n'  # count features
    #Output detailed read assignment results for each read (or fragment if paired end). They are saved to a tab-delimited
    # file that contains four columns including read name, status(assigned or the reason if not assigned), name of target
    # feature/metafeature qqand total number of hits if the read/fragment is counted multiple times.
    #Filtered features: Unassigned_Ambiguity, Unassigned_MultiMapping, Unassigned_NoFeatures, Unassigned_Unmapped

    #combine multimate with unmate
    #alignment

    filename =name+'.sh'
    if iternum == 0:
        newfile = open(filename, 'w')
        newfile.write(
            '#!/bin/bash \n #SBATCH --job-name=serial \n #SBATCH --partition=super \n #SBATCH --nodes=1 \n #SBATCH --time=10-00:00:00 \n #SBATCH --output=./sbatch_output_%j\n')
    if iternum < 40:
        newfile.write(cmd)
        iternum = iternum + 1
        continue
    newfile.close()
    iternum = 0
index.close()
#SAM9267229_1_R1.fastq.gzAligned.out.sam test

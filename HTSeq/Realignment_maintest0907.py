#! usr/bin/env python
import os
import sys
import re
Align_index=sys.argv[1]
Align_HERVRef= sys.argv[2]
Align_Ref=sys.argv[3]
Data_path=sys.argv[4]
name_R1=sys.argv[5]
name_R2=sys.argv[6]
Fasta_path=sys.argv[7]
output_prefix=sys.argv[8]
#bamfile received from cmd_featureCOunts
#hard/soft clipping
#bedtools bamfile2fastq
if name_R1.find('/')>=0:
	name_R1=re.split('/',name_R1)[-1]
if name_R2.find('/')>=0:
	name_R2=re.split('/',name_R2)[-1]
	
pathR1=os.path.join(Fasta_path,name_R1)
pathR2=os.path.join(Fasta_path,name_R2)
output_path=os.path.join(Data_path,output_prefix)
out_NSAM=os.path.join(output_path,(output_prefix+'Aligned.out.sam'))
out_Nfeature=os.path.join(output_path,(output_prefix+'Features.txt'))
feature_details=out_NSAM+'.featureCounts'
HSAM=os.path.join(output_path,'unmated_realignmentAligned.out.sam')
CtrlNorm=os.path.join(output_path,(output_prefix+'Features.txt'))
#if not os.path.exists(output_path):
	#os.makedirs(output_path)
#cmd_STARNORM='STAR --runMode alignReads --runThreadN 40 --genomeDir '+Align_index+' --readFilesIn '+pathR1+' '+pathR2+\
        #' --readFilesCommand zcat --outFilterMultimapNmax 200 --outSAMunmapped Within --outSAMorder PairedKeepInputOrder --outReadsUnmapped None --outFileNamePrefix '+\
        #os.path.join(output_path,output_prefix)+'\n'
#os.system(cmd_STARNORM)
#cmd_featureCounts=cmd = 'featureCounts -T 32 -p -f --donotsort -R -M -a ' + Align_Ref + ' -o ' + out_Nfeature + ' ' + out_NSAM + '\n'
#os.system(cmd_featureCounts)
cmd_pyextract=' '.join(['python Realignment_061617.py',Align_HERVRef,Data_path,name_R1,name_R2,Fasta_path,out_NSAM,feature_details,output_prefix,'\n'])
os.system(cmd_pyextract)
cmd_samCount=' '.join(['python Realignment_samfeatureCount.py',HSAM,CtrlNorm])
os.system(cmd_samCount)


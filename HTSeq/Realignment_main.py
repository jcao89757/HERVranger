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
if len(sys.argv)>9:
	test_status=sys.argv[9]
else:
	test_status='serious'
	print('Will remove all intermediate files.')
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
unmated1 = os.path.join(Data_path, output_prefix, (output_prefix + 'Unmapped.out.mate1'))
unmated2 = os.path.join(Data_path, output_prefix, (output_prefix + 'Unmapped.out.mate2'))
if not os.path.exists(output_path):
	os.makedirs(output_path)
	#--limitOutSAMoneReadBytes 1000000 may not be necessary in other systems
cmd_STARNORM='STAR --runMode alignReads --runThreadN 40 --genomeDir '+Align_index+' --readFilesIn '+pathR1+' '+pathR2+\
        ' --readFilesCommand zcat --limitOutSAMoneReadBytes 1000000 --outFilterMultimapNmax 200 --outSAMunmapped Within --outSAMorder PairedKeepInputOrder --outReadsUnmapped None --outFileNamePrefix '+\
        os.path.join(output_path,output_prefix)+'\n'
#os.system(cmd_STARNORM)
print('1st Alignment completed: '+output_prefix)
cmd_featureCounts=cmd = 'featureCounts -T 32 -p -f --donotsort -R -M -a ' + Align_Ref + ' -o ' + out_Nfeature + ' ' + out_NSAM + '\n'
#os.system(cmd_featureCounts)
print('1st featureCounts completed: '+output_prefix)
cmd_pyextract=' '.join(['python Realignment_061617.py',Align_HERVRef,Data_path,name_R1,name_R2,Fasta_path,out_NSAM,feature_details,output_prefix,'\n'])
#os.system(cmd_pyextract)
print('Putative HERV extraction completed: '+output_prefix)
cmd_STAR = ' '.join(['STAR','--runMode','alignReads','--runThreadN','40','--genomeDir',Align_HERVRef,'--readFilesIn',unmated1,unmated2,'--outFilterMultimapNmax','20','--outFilterMismatchNoverLmax','0.05','outReadsUnmapped','None','--outFileNamePrefix',os.path.join(Data_path, output_prefix, 'unmated_realignment')])
#cmd_STAR = 'STAR --runMode alignReads --runThreadN 40 --genomeDir ' + AlignRef + ' --readFilesIn ' + unmated1 + ' ' + unmated2 + \
           #' --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.05 --outReadsUnmapped None --outFileNamePrefix ' + \
           #os.path.join(Data_path, output_prefix, 'unmated_realignment')
#os.system(cmd_STAR)
cmd_samCount=' '.join(['python Realignment_samfeatureCount.py',HSAM,CtrlNorm])
os.system(cmd_samCount)
print('HERV counting completed: '+output_prefix)

#this section is for test
#for i in list([5,20,50]):
#	for j in list([0.05,0.1,0.3]):
#		cmd_STAR = ' '.join(['STAR','--runMode','alignReads','--runThreadN','40','--genomeDir',Align_HERVRef,'--readFilesIn',unmated1,unmated2,'--outFilterMultimapNmax',str(i),'--outFilterMismatchNoverLmax',str(j),'outReadsUnmapped','None','--outFileNamePrefix',os.path.join(Data_path, output_prefix, '_'.join(['unmated_realigntest',str(i),str(j)]))])
#		os.system(cmd_STAR)
#		cmd_samCount=' '.join(['python Realignment_samfeatureCount.py',os.path.join(output_path,('_'.join(['unmated_realigntest',str(i),str(j)])+'Aligned.out.sam')),CtrlNorm])
#		os.system(cmd_samCount)
#test ends at here

if test_status.find('serious')>=0:
	cmd_rm=' '.join(['rm',out_NSAM,out_Nfeature,feature_details,HSAM,unmated1,unmated2])
	os.system(cmd_rm)
elif test_status.find('test')>=0:
	pass
else:
	print('Wrong test status! Would not remove intermediate files.')
print('You have been terminated: '+output_prefix)


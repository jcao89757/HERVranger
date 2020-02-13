#! usr/bin/env python
#This file is a supplement for Realignment_061617, to provide functions for sequence unmapped extraction.
# Could try extract file from fastq line
# finilized in July 13, 2017, tried in test_combinedict.py
import re
import os
from itertools import izip_longest
class fa_feature(object):

    def __init__(self,patStat,faStat):
        self.patStat=patStat
        self.faStat=faStat

    def isfilled(self):
        if self.patStat is None or self.faStat is None:
            return 0
        else:
            if self.patStat==0:
                return 1
            else:
                return 2

def iter_fasta(iter_items,dict_unsolved,unmated_edit1,unmated_edit2):
    j=0
    for (line1, line2, line3, line4) in iter_items:
        j+=1
        if j%100==0:
            break
        else:
            tmp_fa = re.split('[\n|\t| |@]', line1[0])
            if not dict_unsolved.has_key(tmp_fa[1]):
                dict_unsolved[tmp_fa[1]]=fa_feature(None,None)
            dict_unsolved[tmp_fa[1]].faStat=[line1,line2,line3,line4]
            if dict_unsolved[tmp_fa[1]].isfilled() != 0:
                if dict_unsolved[tmp_fa[1]].isfilled() == 1:
                    write_fasta(dict_unsolved[tmp_fa[1]].faStat,unmated_edit1, unmated_edit2)
                del dict_unsolved[tmp_fa[1]]
    return dict_unsolved

def write_fasta(filled_fa_fea,unmated_edit1,unmated_edit2):
    [line1,line2,line3,line4]=filled_fa_fea
    unmated_edit1.write(line1[0])
    unmated_edit1.write(line2[0])
    unmated_edit1.write(line3[0])
    unmated_edit1.write(line4[0])
    unmated_edit2.write(line1[1])
    unmated_edit2.write(line2[1])
    unmated_edit2.write(line3[1])
    unmated_edit2.write(line4[1])

def seqExtract(Fasta_path, Data_path,name_R1,name_R2, feature_details, unmated1, unmated2,output_prefix,mode):
    pattern = 'initial'
    dict_feature = {'initial': 1, }
    dict_unsolved = {}
    #Copy both fasta files from original path
    Fasta_R1gz=os.path.join(Fasta_path,name_R1)
    #name_tmp=re.split('[_|.]',name)
    #Fasta_R2name='.'.join(['_'.join([name_tmp[0],name_tmp[1],'R2']),'fastq','gz'])
    Fasta_R2gz=os.path.join(Fasta_path,name_R2)
    current_folder=os.path.join(Data_path,output_prefix)
    cmd_cpfastaR1=' '.join(['cp',Fasta_R1gz,current_folder])
    cmd_cpfastaR2=' '.join(['cp',Fasta_R2gz,current_folder])
    Fasta_R1_cur_gz=os.path.join(current_folder,name_R1)
    Fasta_R2_cur_gz=os.path.join(current_folder,name_R2)
    if not os.path.isfile(Fasta_R1_cur_gz):
        os.system(cmd_cpfastaR1)
    cmd_unzip1=' '.join(['gzip','-d',Fasta_R1_cur_gz])
    os.system(cmd_unzip1)
    Fasta_R1_cur_gz=os.path.join(current_folder,name_R1.replace('.gz',''))
    if mode=='paired':
        if not os.path.isfile(Fasta_R2_cur_gz):
            os.system(cmd_cpfastaR2)
        cmd_unzip2=' '.join(['gzip','-d',Fasta_R2_cur_gz])
        os.system(cmd_unzip2)
        Fasta_R2_cur_gz = os.path.join(current_folder, name_R2.replace('.gz', ''))
    else:
        Fasta_R2_cur_gz = os.path.join(current_folder, name_R2.replace('.gz', ''))
        cmd_cpunzip=' '.join(['cp',Fasta_R1_cur_gz,Fasta_R2_cur_gz])
        os.system(cmd_cpunzip)
    unmated_edit1 = open(unmated1, 'w')
    unmated_edit2 = open(unmated2, 'w')
    # Iteration: feature_index
    feature_index = open(feature_details, 'r') 
    fa_R1 = open(Fasta_R1_cur_gz, 'r')
    fa_R2 = open(Fasta_R2_cur_gz, 'r')
    iter_items = izip_longest(*[izip_longest(fa_R1, fa_R2)] * 4)
    i = 0

    for line in feature_index:
        i += 1
        if i % 100 == 0:
            dict_unsolved = iter_fasta(iter_items, dict_unsolved,unmated_edit1,unmated_edit2)

        tmp = re.split('[\n|\t]', line)
        cur = tmp[0]
        if not dict_feature.has_key(cur):  ## deal with the last feature
            if not dict_unsolved.has_key(pattern):
                dict_unsolved[pattern] = fa_feature(None, None)
            dict_unsolved[pattern].patStat = dict_feature[pattern]
            if dict_unsolved[pattern].isfilled() != 0:
                if dict_unsolved[pattern].isfilled() == 1:
                    write_fasta(dict_unsolved[pattern].faStat, unmated_edit1, unmated_edit2)
                del dict_unsolved[pattern]
            dict_feature = {}  # remove last feature
        if re.match(tmp[1], 'Unassigned_NoFeatures') or re.match(tmp[1], 'Unassigned_Unmapped'):
            if dict_feature.has_key(tmp[0]):
                continue
            else:
                pattern = tmp[0]
                dict_feature[tmp[0]] = 0
        else:
            pattern = tmp[0]
            dict_feature[tmp[0]] = 1
            continue

            # deal with the rest of seqs
    for (line1, line2, line3, line4) in iter_items:
        tmp_fa = re.split('[\n|\t| |@]', line1[0])
        if not dict_unsolved[tmp_fa[1]]:
            dict_unsolved[tmp_fa[1]] = fa_feature(None, None)
        dict_unsolved[tmp_fa[1]].faStat = [line1, line2, line3, line4]
        if dict_unsolved[tmp_fa[1]].isfilled() != 0:
            if dict_unsolved[tmp_fa[1]].isfilled() == 1:
                write_fasta(dict_unsolved[tmp_fa[1]].faStat,unmated_edit1, unmated_edit2)
            del dict_unsolved[tmp_fa[1]]

    unmated_edit1.close()
    unmated_edit2.close()
    # Required in this case,to remove copied fasta files

    cmd_rmcp=' '.join(['rm',Fasta_R1_cur_gz,Fasta_R2_cur_gz])
    os.system(cmd_rmcp)

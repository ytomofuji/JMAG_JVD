#!/usr/bin/env python

import sys
import numpy as np
import re
args = sys.argv
INPUT=args[1]
REF=args[2]
FASTA_OUT=args[3]

fasta=[]
fr=open(REF)
lines=fr.readlines()
ID_list=[]
NEWID_list=[]
ID=""
for i in range(1,len(lines)):
    tmp=lines[i].rstrip('\n').split("\t")
    ID_list.append(tmp[0])
    NEWID_list.append(tmp[1])

ref=dict(zip(ID_list,NEWID_list))
f=open(INPUT)
lines=f.readlines()
#prophages
NUM=int(len(lines)/2)
for i in range(NUM):
    NAME=lines[i*2].rstrip('\n')
    SEQ=lines[i*2+1].rstrip('\n')
    tmp=NAME.lstrip('>')
    NEW_NAME=ref[tmp]
    #processed fasta
    fasta.append(">"+NEW_NAME)
    fasta.append(SEQ)

with open(FASTA_OUT,mode='w') as f:
    f.writelines('\n'.join(fasta)+'\n')
#!/usr/bin/env python3
### 2022.09.27



import pysam
from optparse import OptionParser
parser=OptionParser()
parser.add_option('-o', dest='input1')
parser.add_option('-r', dest='input2')
(options, args) = parser.parse_args()
input1=options.input1
input2=options.input2
inputFile1 = pysam.AlignmentFile(input1,"r")
inputFile2 = pysam.AlignmentFile(input2,"r")

def complseq(lt1):
  lt2 = lt1[::-1]
  lt3=''
  for x in lt2:
    if x == 'A': lt3 = lt3 + 'T'
    elif x == 'T': lt3 = lt3 + 'A'
    elif x == 'C': lt3 = lt3 + 'G'
    elif x == 'G': lt3 = lt3 + 'C'
    else: lt3 = lt3 + x
  return lt3

ori_dict = dict()
length_dict = dict()
tuple_dict = dict()
output_dict = dict()
ALU_dict = dict()
for read in inputFile1:
    if read.reference_name == "ALU":
        ALU_dict[read.query_name] = read.to_string()
    else:
        ori_dict[read.query_name] = read.to_string()
        length_dict[read.query_name] = read.query_length
        tuple_dict[read.query_name] = read.cigartuples
        output_dict[read.query_name] = read.to_string()
for read in inputFile2:
    oriRead = ori_dict[read.query_name].split("\t")
    reAlignRead = read.to_string().split("\t")
    oriRead_flag = int(oriRead[1])
    reAlignRead_flag = int(reAlignRead[1])
    if reAlignRead_flag == 4: # realignment not mapped
        continue
    elif (oriRead_flag & 16 == 16 and reAlignRead_flag & 16 == 16) or (oriRead_flag & 16 != 16 and reAlignRead_flag & 16 != 16):# origin alignment and realignment are mapped to the same strand
        modified_Slength = length_dict[read.query_name] - sum([read.cigartuples[i][1] for i in range(len(read.cigartuples)) if read.cigartuples[i][0] in [0,1,4]])
        if tuple_dict[read.query_name][0][0] == 4:
            modified_cigar = read.cigarstring + str(modified_Slength) + "S"
        else:
            modified_cigar = str(modified_Slength) + "S" + read.cigarstring
        reAlignRead[5] = modified_cigar
        reAlignRead[9] = oriRead[9]
        reAlignRead[10] = oriRead[10]
        reAlignRead.append("RA:Z:Y")
        output_dict[read.query_name] = "\t".join(reAlignRead)

    else: # origin alignment and realignment are mapped to different strands 
        modified_Slength = length_dict[read.query_name] - sum([read.cigartuples[i][1] for i in range(len(read.cigartuples)) if read.cigartuples[i][0] in [0,1,4]])
        if tuple_dict[read.query_name][0][0] == 4:
            modified_cigar = str(modified_Slength) + "S" + read.cigarstring
        else:
            modified_cigar = read.cigarstring + str(modified_Slength) + "S"
        reAlignRead[5] = modified_cigar
        reAlignRead[9] = complseq(oriRead[9])
        reAlignRead[10] = oriRead[10][::-1]
        reAlignRead.append("RA:Z:Y")
        output_dict[read.query_name] = "\t".join(reAlignRead)
for read in sorted(output_dict.keys()):
    print(ALU_dict[read])      
    print(output_dict[read])
    
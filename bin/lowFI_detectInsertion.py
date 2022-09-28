#!/usr/bin/env python3
### 2022.09.27
### 1. 


import pysam


inputFile = pysam.AlignmentFile("-","r")

NR=1
marker1 = 0
marker2 = 0

for read in inputFile:
    if(NR%4==1):
        read1TOrd = [read.cigartuples[i][0] for i in range(len(read.cigartuples)) if read.cigartuples[i][0] != 1 and read.cigartuples[i][0] != 2]
        read1TFlag = read.flag
    if(NR%4==2):
        IorDIndexList = [i for i in range(len(read.cigartuples)) if read.cigartuples[i][0] == 1 or read.cigartuples[i][0] == 2]
        if(sum([read.cigartuples[i][1] for i in IorDIndexList]) <= 10):
            read1Line = read.tostring()
            read1Chr = read.reference_name
            read1GStart = read.reference_start # sam, 1-based
            read1GEnd = read.reference_end
            read1GOrd = [read.cigartuples[i][0] for i in range(len(read.cigartuples)) if read.cigartuples[i][0] != 1 and read.cigartuples[i][0] != 2]
            read1GM = max([read.cigartuples[i][1] for i in range(len(read.cigartuples)) if read.cigartuples[i][0] == 0 ])
            read1GFlag = read.flag
            marker1 = 1
        else:
            marker1 = 0
    if(NR%4==3):
        read2TOrd = [read.cigartuples[i][0] for i in range(len(read.cigartuples)) if read.cigartuples[i][0] != 1 and read.cigartuples[i][0] != 2]
        read2TFlag = read.flag
    if(NR%4==0):
        IorDIndexList = [i for i in range(len(read.cigartuples)) if read.cigartuples[i][0] == 1 or read.cigartuples[i][0] == 2]
        if(sum([read.cigartuples[i][1] for i in IorDIndexList]) <= 10):
            read2Line = read.tostring()
            read2Chr = read.reference_name
            read2GStart = read.reference_start # sam, 1-based
            read2GEnd = read.reference_end
            read2GOrd = [read.cigartuples[i][0] for i in range(len(read.cigartuples)) if read.cigartuples[i][0] != 1 and read.cigartuples[i][0] != 2]
            read2GM = max([read.cigartuples[i][1] for i in range(len(read.cigartuples)) if read.cigartuples[i][0] == 0 ])
            read2GFlag = read.flag
            marker2 = 1
        else:
            marker2 = 0
    if(NR%4==0 and marker1 == 1 and marker2 ==1 and read1Chr == read2Chr):
        if(read2GStart >= read1GStart):
            if(read1GFlag & 16 != 16 and read2GFlag & 16 == 16 and read1TFlag == 0 and read2TFlag == 16 and read1GOrd[0] == 0 and read1GOrd[-1] == 4 and read2GOrd[0] == 4 and read2GOrd[-1] == 0 and read1TOrd[0] == 4 and read1TOrd[-1] == 0 and read2TOrd[0] == 0 and read2TOrd[-1] == 4): # ALU + strand 
                if(abs(read1GEnd - read2GStart) <= 50):
                    print(read2Chr + "\t" + str(read2GStart) + "\t" + str(read2GStart) + "\t" + read2Line.split("\t")[0].split("_")[0])
            if(read1GFlag & 16 != 16 and read2GFlag & 16 == 16 and read1TFlag == 16 and read2TFlag == 0 and read1GOrd[0] == 0 and read1GOrd[-1] == 4 and read2GOrd[0] == 4 and read2GOrd[-1] == 0 and read1TOrd[0] == 0 and read1TOrd[-1] == 4 and read2TOrd[0] == 4 and read2TOrd[-1] == 0): # ALU - strand
                if(abs(read1GEnd - read2GStart) <= 50):
                    print(read2Chr + "\t" + str(read2GStart) + "\t" + str(read2GStart) + "\t" + read2Line.split("\t")[0].split("_")[0])
        else:
            if(read2GFlag & 16 != 16 and read1GFlag & 16 == 16 and read2TFlag == 0 and read1TFlag == 16 and read2GOrd[0] == 0 and read2GOrd[-1] == 4 and read1GOrd[0] == 4 and read1GOrd[-1] == 0 and read2TOrd[0] == 4 and read2TOrd[-1] == 0 and read1TOrd[0] == 0 and read1TOrd[-1] == 4): # ALU + strand 
                if(abs(read2GEnd - read1GStart) <= 50):
                    print(read1Chr + "\t" + str(read1GStart) + "\t" + str(read1GStart) + "\t" + read1Line.split("\t")[0].split("_")[0])
            if(read2GFlag & 16 != 16 and read1GFlag & 16 == 16 and read2TFlag == 16 and read1TFlag == 0 and read2GOrd[0] == 0 and read2GOrd[-1] == 4 and read1GOrd[0] == 4 and read1GOrd[-1] == 0 and read2TOrd[0] == 0 and read2TOrd[-1] == 4 and read1TOrd[0] == 4 and read1TOrd[-1] == 0): # ALU - strand
                if(abs(read2GEnd - read1GStart) <= 50):
                    print(read1Chr + "\t" + str(read1GStart) + "\t" + str(read1GStart) + "\t" + read1Line.split("\t")[0].split("_")[0])               
    NR = NR + 1
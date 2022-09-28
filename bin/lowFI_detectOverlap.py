#!/usr/bin/env python3
### detectDifLoc.py
### 2022.09.27



import pysam


class ReadGT:
    """
    Class for calculation of the overlapping length and proportion it constitutes in TE consensus sequence mapped part. High proportion (i.e., no less than 50%) indicates the necessity of realignment of the soft-clipped part. The length of soft-clipped part which overlaps with ALU mapped part should be no more than 50%.
    """
    def __init__(self,Gcigartuple,Tcigartuple):
        self.Gcigartuple = Gcigartuple
        self.Tcigartuple = Tcigartuple
        self.GTcigarlist = self.genList()
        self.GTMatchOverlap = len([i for i in self.GTcigarlist if i == ["0","0"]])/len([i for i in self.GTcigarlist if i[1] == "0"])
        # self.GTSoftOverlap = len([i for i in self.GTcigarlist if i == ["4","4"]])/len([i for i in self.GTcigarlist if i[0] == "4"])
    def genList(self):
        Gcigarstring = "".join([str(self.Gcigartuple[i][0]) * int(self.Gcigartuple[i][1]) for i in range(len(self.Gcigartuple)) if self.Gcigartuple[i][0] in [0,1,4]])
        Tcigarstring = "".join([str(self.Tcigartuple[i][0]) * int(self.Tcigartuple[i][1]) for i in range(len(self.Tcigartuple)) if self.Tcigartuple[i][0] in [0,1,4]])
        return [[i,j] for i,j in zip(Gcigarstring, Tcigarstring)]

inputFile = pysam.AlignmentFile("-","r")
NR = 1
for read in inputFile:
    if(NR%2==1):
        readTcigartuple=read.cigartuples
    else:
        readGcigartuple=read.cigartuples
        readGT = ReadGT(readGcigartuple,readTcigartuple)
        if(readGT.GTMatchOverlap >= 0.5):
            print(read.to_string())
    NR = NR + 1

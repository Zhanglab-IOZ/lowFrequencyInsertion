#!/usr/bin/env python3
###### convertFastq.py ######
### 2022.09.27
### 1. Convert read - strand to its reverse complementary sequence.

from optparse import OptionParser
parser=OptionParser()
parser.add_option('-i', dest='inputFastq')
(options, args) = parser.parse_args()
fi=options.inputFastq

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

if __name__ == '__main__':
  fin=open(fi,'r')
  lin=fin.readlines()
  for i in range(len(lin)):
    if (i+1)%4 == 1:
        print(lin[i],end="")
        flag = int(lin[i].split("_")[1])
    if (i+1)%4 == 3:
        print("+")
    if (i+1)%4 == 2:
        if flag & 0x10 == 16:
            print(complseq(lin[i].strip())+"\n",end="")
        else:
            print(lin[i],end="")
    if (i+1)%4 == 0:
        if flag & 0x10 == 16:
            print(lin[i].strip()[::-1]+"\n",end="")
        else:
            print(lin[i],end="")
  fin.close()

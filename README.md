# lowFrequencyInsertion  
## Dependency  
novoalign version: 3.09.04  
parallel<sup>[[1]](https://github.com/MarcelloMalpighi/lowFrequencyInsertion/edit/main/README.md#references)</sup> version: 20220722  
samtools version: 1.15.1  
  
## getSoftClipped.sh  
Get primary alignments contains only 1 soft-clipped part and  
the length of this part must be not more/less than input limits.  
Output is composed of a sam file of chosen alignments and  
a corresponging fastq file of soft-clipped parts.  
Usage:  
[-i <input file, bam/sam>]  
[-o <output file name, suffix will be added automatically>]  
[-u <upper limit of soft-clipped part length, limit itself is included>]  
[-l <lower limit of soft-clipped part length, limit itself is included>]  
[-h \<help>]  
  
## reAlignRef.sh  
Realign soft-clipped parts to TE consensus sequences with novoalign.  
Output is composed a sorted.bam file containing all mapped soft-clipped  
parts sorted by coordinates and a novoalign report. Parallel execution is allowed  
with GNU Parallel.  
Usage:  
[-i <input file, fastq>]  
[-o <output file name, suffix will be added automatically>]  
[-c \<TE consensus sequences novoalign index file>]  
[-r <upper limit of secondary alignments novoalign will report, default: 10>]  
[-R <standard of unique alignments novoalign will use, higher, more stringent, less unique aliments, default: 5>]  
[-p <number of jobs to be run in parallel, default: 2>]  
[-h \<help>]  
  
## References  
[1] Tange, Ole. (2018). GNU Parallel 2018. In GNU Parallel 2018 (p. 112). Ole Tange. https://doi.org/10.5281/zenodo.1146014

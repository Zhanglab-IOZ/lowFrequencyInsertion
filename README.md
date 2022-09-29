# lowFrequencyInsertion  
## Dependencies  
1. novoalign version: 3.09.04  
2. parallel<sup>[[1]](https://github.com/MarcelloMalpighi/lowFrequencyInsertion/edit/main/README.md#references)</sup> version: 20220722  
3. pysam version: 0.19.1  
4. python version: 3.10.5  
5. samtools version: 1.15.1  
  
## Run lowFI
```
lowFI  
Detect low-frequency ALU insertions supported by single soft-clipped read pair.
  
Usage: lowFI [options]  
[-i <input file, the absolute path is necessary, bam/sam, mandatory>]  
[-o <output file name, suffix will be added automatically, mandatory>]  
[-u <upper limit of soft-clipped part length, limit itself is included, optional, default: 130>]  
[-l <lower limit of soft-clipped part length, limit itself is included, optional, default: 20>]  
[-p <number of jobs to be run in parallel, optional, default: 2>]  
[-m <memory per thread used for samtools sort, optional, defalut: 2G>]  
[-T <ALU consensus sequences novoalign index file, mandatory>]  
[-G <Genome novoalign index file, mandatory>]  
[-R <ALU annotation file, bed, mandatory>]  
[-X <nonreference insertion detection result, bed, optional>]  
[-h <help>]  
```
  
## References  
[1] Tange, Ole. (2018). GNU Parallel 2018. In GNU Parallel 2018 (p. 112). Ole Tange. https://doi.org/10.5281/zenodo.1146014

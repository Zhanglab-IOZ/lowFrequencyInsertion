#!/bin/bash
### 2022.08.31
### 1. Filter alignments with flag 4(unmapped) 256(secondary) 1024(duplicates) 2048(supplementary).
### 2. Get alignments contains only 1 soft-clipped part whose length must be not more/less than input limits and their mate alignments.



# usage for output
usage(){ 
echo -e "getSoftClipped\n
Get primary alignments contains only 1 soft-clipped part and the length 
of this part must be not more/less than input limits. Output is composed
of a sam file of chosen alignments with their mate alignments and a 
corresponging fastq file of chosen reads.\n
Usage: 
[-i <input file, bam/sam, mandatory>]
[-o <output file name, suffix will be added automatically, mandatory>]
[-u <upper limit of soft-clipped part length, limit itself is included, optional, default: 130>]
[-l <lower limit of soft-clipped part length, limit itself is included, optional, default: 20>]
[-t <threads used for samtools sort, optional, default: 2>]
[-m <memory per thread used for samtools sort, optional, defalut: 2G>]
[-h <help>]"; 
exit;
}

# get opts
while getopts ":hi:o:u:l:t:m:" opt
do
	case "${opt}" in
	h)usage;;
	i)input=${OPTARG};;
	o)output=${OPTARG};;
	u)upper=${OPTARG};;
	l)lower=${OPTARG};;
    t)threads=${OPTARG};;
    m)memory=${OPTARG};;
	*)
        	echo "Error! Unexpected parameters: $*!" 1>&2
		usage
		;;
	esac
done

[ -z $upper ] && upper=130
[ -z $lower ] && lower=20
[ -z $threads ] && threads=2
[ -z $memory ] && memory="2G"

# prepare sam file
if test $input && test $output && test $upper && test $lower
then
samtools view -u -F 2304 ${input} | samtools sort -n -@ $threads -m $memory -O sam | \
awk -v up=$upper -v low=$lower -f getSoftClipped.awk.sh > ${output}.sam
else
echo "Error! Lack of mandatory parameters!"
usage
fi

# output fastq file with QNAME FLAG CIGAR RNAME POS and sam file
awk '{printf "%s%s%s%s%s%s%s%s%s%s\n","@",$1,"_",$2,"_",$6,"_",$3,"_",$4;printf "%s\n%s\n%s\n",$10,"+",$11}' ${output}.sam > ${output}.fastq
samtools view -H ${input} -o ${output}.header
#echo -e "@PG\tID:lowFrequencyInsertion\tVN:0.0.1\tCL:getSoftClipped" >> ${output}.header
cat ${output}.sam >> ${output}.header
mv ${output}.header ${output}.sam
echo "getSoftClipped done"

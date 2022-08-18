#!/bin/bash
### 2022.08.18
### 1. Filter alignments with flag 4(unmapped) 256(secondary) 1024(duplicates) 2048(supplementary).
### 2. Get alignments contains only 1 soft-clipped part and the length of this part must be not less than 10 bases.



# usage for output
usage(){ 
echo -e "getSoftClipped\n
Get primary alignments contains only 1 soft-clipped part and 
the length of this part must be not more/less than input limits.
Output is composed of a sam file of chosen alignments and
a corresponging fastq file of soft-clipped parts.\n
Usage: 
[-i <input file, bam/sam>]
[-o <output file name, suffix will be added automatically>]
[-u <upper limit of soft-clipped part length, limit itself is included>]
[-l <lower limit of soft-clipped part length, limit itself is included>]
[-h <help>]"; 
exit;
}
# Parallel execution is allowed with GNU Parallel.
# [-p <number of jobs to be run in parallel, default: 2>]

# get opts
while getopts ":hi:o:u:l:" opt
do
	case "${opt}" in
	h)usage;;
	i)input=${OPTARG};;
	o)output=${OPTARG};;
	u)upper=${OPTARG};;
	l)lower=${OPTARG};;	
	*)
        	echo "Error! Unexpected parameters: $*!" 1>&2
		usage
		;;
	esac
done
#    p)parallel=${OPTARG};;
# prepare sam file in parallel
# [ -z $parallel ] && parallel=2

# prepare sam file
if test $input && test $output && test $upper && test $lower
then
regular=(`seq -s "S;" $lower $upper`"S")
ul=$(($upper-$lower+1))
samtools view -F 3332 ${input} | \
awk '{split($6,a,"S");if(length(a)==2){if(a[2]==""){printf "%s\t%s\n",$0,"TL:Z:"}else{printf "%s\t%s\n",$0,"HD:Z:"}}}' | \
awk -v ul=$ul -v re=$regular 'BEGIN{split(re,re_ar,";")}{for(i=1;i<=ul;++i){if(match($6,"^"""re_ar[i]"""|"re_ar[i]"""$")){printf "%s%s\n",$0,re_ar[i];break}}}' > ${output}.sam

### failed parallel with awk errors
# samtools view -F 3332 $input -o ${output}_noheader
# split -d -a 4 -n l/$parallel --additional-suffix=.inter ${output}_noheader ${output}_noheader_
# rm ${output}_noheader
# ls ${output}_noheader_*.inter | parallel --will-cite --linebuffer --progress -j $parallel "mkdir {/}_dir;cd {/}_dir;awk '{split(\$6,a,\"S\");if(length(a)==2){if(a[2]==\"\"){printf \"%s\t%s\n\",\$0,\"TL:Z:\"}else{printf \"%s\t%s\n\",\$0,\"HD:Z:\"}}}' | awk -v ul=$ul -v re=$regular 'BEGIN{split(re,re_ar,\";\")}{for(i=1;i<=ul;++i){if(match(\$6,\"^\"\"\"re_ar[i]\"\"\"|\"re_ar[i]\"\"\"\$\")){printf \"%s%s\n\",\$0,re_ar[i];break}}}' > {/}.sam"
# for i in `ls ${output}_noheader_*.inter`
# do
# rm $i
# done

### an alternative algorithm but slower
# awk -v ul=$ul -v re=$regular 'BEGIN{split(re,re_ar,";")}{split($6,a,"S");if(length(a)==2){if(a[2]==""){for(i=1;i<=ul;++i){if(match($6,re_ar[i]"""$")){printf "%s\t%s%s\n",$0,"TL:Z:",re_ar[i]}}}else{for(i=1;i<=ul;++i){if(match($6,"^"""re_ar[i])){printf "%s\t%s%s\n",$0,"HD:Z:",re_ar[i]}}}}}' > ${output}.sam

else
echo "Error! Lack of mandatory parameters!"
usage
fi

# output fastq file with QNAME FLAG CIGAR RNAME POS and sam file
awk '{printf "%s%s%s%s%s%s%s%s%s%s\n","@",$1,"_",$2,"_",$6,"_",$3,"_",$4;if($NF~/HD/){len = substr($NF,6,length($NF)-6);printf "%s\n%s\n%s\n",substr($10,1,len),"+",substr($11,1,len)}else{len = substr($NF,6,length($NF)-6);printf "%s\n%s\n%s\n",substr($10,length($10)-len+1,len),"+",substr($11,length($11)-len+1,len)}}' ${output}.sam > ${output}.fastq
samtools view -H ${input} -o ${output}.header
#echo -e "@PG\tID:lowFrequencyInsertion\tVN:0.0.1\tCL:getSoftClipped" >> ${output}.header
cat ${output}.sam >> ${output}.header
mv ${output}.header ${output}.sam
echo "getSoftClipped done"
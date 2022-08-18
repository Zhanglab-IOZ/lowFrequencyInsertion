#!/bin/bash
### 2022.08.18
### 1. Realign fastq file of soft-clipped parts with novoalign to TE consensus sequences.
### 2. Allow multiple novoalign at the same time.



# usage for output
usage(){
echo -e "reAlignRef\n
Realign soft-clipped parts to TE consensus sequences with novoalign.
Output is composed a sorted.bam file containing all mapped soft-clipped
parts sorted by coordinates and a novoalign report. Parallel execution is allowed
with GNU Parallel.\n
Usage: 
[-i <input file, fastq>]
[-o <output file name, suffix will be added automatically>]
[-c <TE consensus sequences novoalign index file>]
[-r <upper limit of secondary alignments novoalign will report, default: 10>]
[-R <standard of unique alignments novoalign will use, higher, more stringent, less unique aliments, default: 5>]
[-p <number of jobs to be run in parallel, default: 2>]
[-h <help>]";
exit;
}

# get opts
while getopts ":hi:o:c:r:R:p:" opt
do
    case "${opt}" in
    h)usage;;
    i)input=${OPTARG};;
    o)output=${OPTARG};;
    c)index=${OPTARG};;
    r)upper=${OPTARG};;
    R)unique=${OPTARG};;
	p)parallel=${OPTARG};;
    *)
        echo "Error! Unexpected parameters: $*!" 1>&2
        usage
        ;;
    esac
done

# parallel novoalign
[ -z $upper ] && upper=10
[ -z $unique ] && unique=5
[ -z $parallel ] && parallel=2

if test $input && test $output && test $index
then
total_len=`wc -l ${input} | awk '{print $1}'`
total_len=$(($total_len/4))
batch_len=$((4*($total_len/$parallel+1)))
split -d -a 4 -l $batch_len --additional-suffix=.inter ${input} ${output}_novoalign_
ls ${output}_novoalign_*.inter | parallel --will-cite --progress -j $parallel -k novoalign -d $index -f {} -o SAM -r All $upper -R $unique ">" {}.sam "2>" {}.novoalign.log
else
echo "Error! Lack of mandatory parameters!" 1>&2
usage
fi

# combine, sort and output realignment results
for i in `ls ${output}_novoalign_*.inter`
do
if [ "${i}" == "${output}_novoalign_0000.inter" ]
then
samtools view -h ${i}.sam -o ${output}.sam
cat ${i}.novoalign.log > ${output}.novoalign.log
else
samtools view ${i}.sam >> ${output}.sam
cat ${i}.novoalign.log >> ${output}.novoalign.log
fi
rm ${i}
rm ${i}.sam
rm ${i}.novoalign.log
done
samtools view -u -F 4 ${output}.sam | samtools sort -@ $parallel -m 2G -O bam -o ${output}.sorted.bam # unmapped reads are abandoned
rm ${output}.sam
echo "reAlignRef done"
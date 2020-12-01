#!/bin/bash

inputfile1="/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq"
inputfile2="/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"




if [ ! -f $inputfile1 ]
then
	echo 'ERROR:' $inputfile1 'not found.'
fi
if [ ! -f $inputfile2 ]
then
	echo 'ERROR:' $inputfile2 'not found.'
fi



echo "Input files:" $inputfile1 $inputfile2

extension='.'$(echo $inputfile1 | rev | cut -d. -f1 | rev)

outputfile_sample1=${inputfile1%$extension*}'_sample1_interleaved.fq'
outputfile_sample2=${inputfile1%$extension*}'_sample2_interleaved.fq'



echo "Demultiplexing input files ..."

egrep -h -A 2 -B 1 'GCCACATA|GCGAGTAA' $inputfile1 $inputfile2 | grep -v -- -- > $outputfile_sample1
egrep -h -A 2 -B 1 'GAGCTGAA|GATAGACA' $inputfile1 $inputfile2 | grep -v -- -- > $outputfile_sample2




echo "Sorting demultiplexed files ..."

cat $outputfile_sample1 | paste - - - - | sort -k1,1 -S 2G | tr '\t' '\n' > ${outputfile_sample1%$extension*}'_sorted.fq'
cat $outputfile_sample2 | paste - - - - | sort -k1,1 -S 2G | tr '\t' '\n' > ${outputfile_sample2%$extension*}'_sorted.fq'

echo "... demulitplexed files sorted."






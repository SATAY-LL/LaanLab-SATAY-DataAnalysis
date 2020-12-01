#!/bin/bash

inputfile1="/home/gregoryvanbeek/Documents/data_processing/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq"
inputfile2="/home/gregoryvanbeek/Documents/data_processing/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"
#inputfile1="/home/gregoryvanbeek/Documents/data_processing/test1.fq"
#inputfile2="/home/gregoryvanbeek/Documents/data_processing/test2.fq"

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

outputfile_sample1_forward=${inputfile1%$extension*}'_sample1_forward_container.fq'
outputfile_sample1_reverse=${inputfile1%$extension*}'_sample1_reverse_container.fq'
outputfile_sample2_forward=${inputfile1%$extension*}'_sample2_forward_container.fq'
outputfile_sample2_reverse=${inputfile1%$extension*}'_sample2_reverse_container.fq'



echo "Demultiplexing input files ..."



#awk '/test/' $inputfile1 > $outputfile_sample1_forward
#sed -n -e '/test/{x;1!p;g;$!N;p;D;}' -e h $inputfile1
#grep -h -A 2 -B 1 'testf' $inputfile1 $inputfile2 | grep -v -- -- > $outputfile_sample1_forward
#grep -h -A 2 -B 1 'testr' $inputfile1 $inputfile2 | grep -v -- -- > $outputfile_sample1_reverse
grep -h -A 2 -B 1 'GCCACATA' $inputfile1 $inputfile2 | grep -v -- -- > $outputfile_sample1_forward
grep -h -A 2 -B 1 'GCGAGTAA' $inputfile1 $inputfile2 | grep -v -- -- > $outputfile_sample1_reverse
grep -h -A 2 -B 1 'GAGCTGAA' $inputfile1 $inputfile2 | grep -v -- -- > $outputfile_sample2_forward
grep -h -A 2 -B 1 'GATAGACA' $inputfile1 $inputfile2 | grep -v -- -- > $outputfile_sample2_reverse



echo "... demultiplexing complete."
echo "Re-pairing output files ..."
bbmappath=/home/gregoryvanbeek/Documents/Software/BBMap/bbmap/

bash ${bbmappath}repair.sh -Xmx1g in=$outputfile_sample1_forward in2=$outputfile_sample1_reverse out=${inputfile1%$extension*}'_sample1_forward.fq' out2=${inputfile1%$extension*}'_sample1_reverse.fq'
rm $outputfile_sample1_forward
rm $outputfile_sample1_reverse
bash ${bbmappath}repair.sh -Xmx1g in=$outputfile_sample2_forward in2=$outputfile_sample2_reverse out=${inputfile1%$extension*}'_sample2_forward.fq' out2=${inputfile1%$extension*}'_sample2_reverse.fq'
rm $outputfile_sample2_forward
rm $outputfile_sample2_reverse

echo "... re-pairing complete."

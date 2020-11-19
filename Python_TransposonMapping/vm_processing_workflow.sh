#!/bin/bash

### This workflow is developed for automatically preprocessing sequencing data.
### The variables stored in the USER SETTINGS block should be checked and potentially changed for each dataset.
### The other lines should be left as they are.
### Within the folder where the datafile is located, three more folders are generated for the output of the quality report, trimming and alignment.
### Settings for trimming and alignment should be set in the 'trimming_settings' and 'alignment_settings' variables respectively.
### Adapter sequences for trimming should be stored in '~/Documents/Software/BBMap/bbmap/resources/adapters.fa' respecting fasta convention.
### The software is called in this order: 1.Quality checking (fastqc) 2.trimming (bbduk or trimmomatic) 3.Quality checking trimmed data (fastqc) 4.alignment (bwa) 5.converting sam to bam (samtools) 6.indexing bam file (sambamba) 7.transposon mapping (python).
### The quality checking (either before or after trimming) can be disabled.
###
### There are two trimming software packages, 'bbduk' and 'trimmomatic'.
### You can use either one of them by settings the appropriate option in the USER SETTINGS.
### When selecting one trimming software, the options of the other program will be ignored.
###
### The workflow will ask after the first quality report if you want to continue (if enabled).
### Pressing 'n' will abort the workflow and allows you to make changes according to the quality report.
### The quality report can be accessed at the location [yourdatafolder]/fastqc_out/[filename].html
### Restarting the workflow cause it to skip over the initial quality report (unless you deleted it) and continues with the trimming and alignment steps.
###
### In case of emergency, press ctrl-c (possibly multiple times) in the terminal.





####################### USER SETTINGS ######################
# Define whether data is paired-end ('t' for paired-end, 'f' for single end)
paired=F


# Define filename (can also be a zipped file ending with .gz). Use filename2 for paired end or leave empty for single end or interleaved paired end (i.e. paired end reads are in single file).
filepath1=/home/gregoryvanbeek/Documents/data_processing/testfolder/SRR062634.filt.fastq.gz
filepath2=''


###### Set options for trimming software ######
# set 'b' for bbduk, set 't' for trimmomatic
trimming_software='b'

###    bbduk    ###
trimming_settings_bbduk='ktrim=r k=4 hdist=0 qtrim=r trimq=4'
## Set adapter sequences
## Open file using xdg-open ~/Documents/Software/BBMap/bbmap/resources/adapters.fa
###################

### trimmomatic ###
trimmomatic_initialization='-phred33'
trimming_settings_trimmomatic='ILLUMINACLIP:adapters.fa:0:30:10 SLIDINGWINDOW:10:4 MINLEN:30'
## Set adapter sequences
## Open file using xdg-open ~/Documents/Software/BBMap/bbmap/resources/adapters.fa
###################
###############################################


# Set options for alignment software (bwa mem)
alignment_settings='-B 2 -O 3'


# Create sorted and indexed bam file ('y' for yes, 'n' for no)?
sort_and_index=T


# Apply transposon mapping (requires sort_and_index=T)
mapping=T


# Save sam file ('y' for yes, 'n' for no)? This file is always converted to its binary equivalent (.bam ) and the sam file is rarely used but takes up relatively a lot of memory.
delete_sam=F


# Open adapters.fa file after the first quality check in order to change the adapters for trimming.
#open_adapters=F


# Create quality report of raw data (before trimming)?
quality_check_raw=F


# Create quality report of trimmed data (after trimming)?
quality_check_trim=F


# Determine whether the script should automatically continue after creating the first quality report. Set to yes if you might want to make changes depending on the quality report of the raw data.
qualitycheck_interrupt=F

############################################################





echo 'Preparing processing for' $(basename ${filepath1}) '...'
echo ''


pathdata=$(dirname ${filepath1})
filename1=$(basename ${filepath1})
if [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filename2} ]]
then
	filename2=$(basename ${filepath2})
fi



if [ ! -f ${filepath1} ]
then
	echo 'ERROR: File' ${filepath1} 'not found.' && exit 1
fi


if [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filename2} ]]
then
	if [ ! -f ${filepath2} ]
	then
		echo 'ERROR: File' ${filepath2} 'not found.' && exit 1
	fi
fi



# Get extension of the file
extension='.'$(echo $filename1 | rev | cut -d. -f1 | rev)
if [[ extension='.gz' ]]
then
	extension='.'$(echo $filename1 | rev | cut -d. -f2 | rev)
fi



# Define filename for trimming and alignment results
filename_trimmed1=${filename1%$extension*}'_trimmed.fastq'
if ! [[ -z ${filename2} ]] #if not filename2 is empty string
then
	filename_trimmed2=${filename2%$extension*}'_trimmed.fastq'
fi
filename_sam=${filename1%$extension*}'_trimmed.sam'
filename_bam=${filename1%$extension*}'_trimmed.bam'
filename_sort=${filename1%$extension*}'_trimmed.sorted.bam'


# Define path output directory fastqc
if [[ ${quality_check_raw} =~ ^[tT]$ ]] || [[ ${quality_check_trim} =~ ^[tT]$ ]]
then
	path_fastqc_out=${pathdata}/fastqc_out
	[ ! -d ${path_fastqc_out} ] && echo 'Creating fastqc output folder ...' && mkdir ${path_fastqc_out} || echo 'Folder for fastqc output exists with name:' $(basename ${path_fastqc_out})
fi

# Define path output directory trimming
path_trimm_out=${pathdata}/trimm_out
[ ! -d ${path_trimm_out} ] && echo 'Creating trimming output folder ...' && mkdir ${path_trimm_out} || echo 'Folder for trimming output exists with name:' $(basename ${path_trimm_out})

# Define path output directory alignment
path_align_out=${pathdata}/align_out
[ ! -d ${path_align_out} ] && echo 'Creating alignment output folder ...' && mkdir ${path_align_out} || echo 'Folder for alignment output exists with name:' $(basename ${path_align_out})

# Define path shared folder
path_sf=/media/sf_VMSharedFolder_Ubuntu64_1/

# Define paths to reference genomes (both S288C and W303)
path_refgenome=/home/gregoryvanbeek/Documents/Reference_Sequences/Reference_Sequence_S288C/S288C_reference_sequence_R64-2-1_20150113.fsa
name_refgenome='S288C'
if [ ! -f ${path_refgenome} ] #if path to reference genome does not exist
then
	echo 'ERROR: Reference genome not found at location:' ${path_refgenome} && exit 1
else
	echo 'Reference genome:' ${name_refgenome}
fi

# Define path bbduk software
path_bbduk_software=~/Documents/Software/BBMap/bbmap/
path_bbduk_adapters=${path_bbduk_software}/resources/adapters.fa
[ ! -d ${path_ddbuk_software} ] && echo 'WARNING: Path to bbduk software does not exists.'

# Define path trimmomatic software
path_trimm_software=~/Documents/Software/Trimmomatic-0.39/
[ ! -d ${path_trimm_software} ] && echo 'WARNING: Path to trimmomatic software does not exists.'

# Define path to python script
path_python_codes=~/Documents/Software/python_codes/
[ ! -d ${path_python_codes} ] && echo 'WARNING: Path to python codes does not exists.'






### Start processing workflow
echo ''
echo 'Start processing ...'
echo ''


# Quality checking raw data
if [[ ${quality_check_raw} =~ ^[tT]$ ]]
then
	if [[ ! -e ${path_fastqc_out}/${filename1%$extension*}'_fastqc.html' ]]
	then
		echo 'Quality checking raw data ...'
		fastqc --outdir ${path_fastqc_out} ${pathdata}/${filename1}
		echo 'Quality checking raw data completed. Results are stored at' ${path_fastqc_out}
		echo ''

		if [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filename2} ]]
		then
			fastqc --outdir ${path_fastqc_out} ${pathdata}/${filename2}
			echo 'Quality checking raw data paired end reads completed. Results are stored at' ${path_fastqc_out}
			echo ''
		fi
	else
		echo 'Quality report raw data already exists. Skipping fastqc'
	fi


	if [[ ${qualitycheck_interrupt} =~ ^[tT]$ ]]
	then
		read -p 'Continue processing? (press "y" if yes, press "n" if no): ' -n 1 -r
		echo
		if [[ ! $REPLY =~ ^[yY]$ ]]
		then
			exit 1
		fi
	fi
fi


#if [[ ${open_adapters} =~ ^[tT]$ ]]
#then
#	echo "Adapter.fa file is being opened..."
#	xdg-open ~/Documents/Software/BBMap/bbmap/resources/adapters.fa
#	read -s -p "Press enter to continue"
#fi


# Trimming
if [[ ${trimming_software} =~ ^[bB]$ ]]
then
	if [[ ${paired} =~ ^[fF]$ ]]
	then
		echo 'Data trimming using bbduk single end reads...'
		${path_bbduk_software}bbduk.sh -Xmx1g in=${pathdata}/${filename1} out=${path_trimm_out}/${filename_trimmed1} ref=${path_bbduk_adapters} ${trimming_settings_bbduk}
		echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1}
		echo ''
	elif [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filename2} ]]
	then
		echo 'Data trimming using bbduk paired end reads...'
		${path_bbduk_software}bbduk.sh -Xmx1g in1=${pathdata}/${filename1} out1=${path_trimm_out}/${filename_trimmed1} in2=${pathdata}/${filename2} out2=${path_trimm_out}/${filename_trimmed2} ref=${path_bbduk_adapters} ${trimming_settings_bbduk}
		echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1} 'and for the paired end reads in' ${path_trimm_out}/${filename_trimmed2}
		echo ''

	elif [[ ${paired} =~ ^[tT]$ ]] && [[ -z ${filename2} ]]
	then
		echo 'Data trimming using bbduk paired end reads...'
		${path_bbduk_software}bbduk.sh -Xmx1g interleaved=t in=${pathdata}/${filename1} out=${path_trimm_out}/${filename_trimmed1} ref=${path_bbduk_adapters} ${trimming_settings_bbduk}
		echo 'Trimming with bbduk is completed. Results are stored in' ${path_trimm_out}/${filename_trimmed1}
		echo ''
	fi

elif [[ ${trimming_software} =~ ^[tT]$ ]]
then
	if [[ ${paired} =~ ^[fF]$ ]]
	then
		echo 'Data trimming using trimmomatic ...'
		currentpath=$(pwd)
		cd ${path_bbduk_software}/resources/
		java -jar ${path_trimm_software}trimmomatic-0.39.jar SE ${trimmomatic_initialization} ${pathdata}/${filename1} ${path_trimm_out}/${filename_trimmed1} ${trimming_settings_trimmomatic}
		cd ${currentpath}

	elif [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filename2} ]]
	then
		echo 'Data trimming using trimmomatic ...'
		currentpath=$(pwd)
		cd ${path_bbduk_software}/resources/
		java -jar ${path_trimm_software}trimmomatic-0.39.jar PE ${trimmomatic_initialization} ${pathdata}/${filename1} ${pathdata}/${filename2} ${path_trimm_out}/${filename_trimmed1} ${path_trimm_out}/${filename_trimmed1%_trimmed.fastq*}'_trimmedorphanedreads.fastq' ${path_trimm_out}/${filename_trimmed2} ${path_trimm_out}/${filename_trimmed1%_trimmed.fastq*}'_trimmedorphanedreads.fastq' ${trimming_settings_trimmomatic}
		cd ${currentpath}

	elif [[ ${paired} =~ ^[tT]$ ]] && [[ -z ${filename2} ]]
	then
		echo 'Enter two input files for using paired end reads with Trimmomatic.'
		exit 1
	fi

else
	echo 'Trimming software not recognized, please check settings' && exit 1
fi


# Quality report trimmed data
if [[ ${quality_check_trim} =~ ^[tT]$ ]]
then
	echo 'Quality checking trimmed data ...'
	fastqc --outdir ${path_fastqc_out} ${path_trimm_out}/${filename_trimmed1}
	echo 'Quality checking trimmed data completed. Results are stored at' ${path_fastqc_out}
	echo ''
	if [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filename2} ]]
	then
		echo 'Quality checking trimmed data paired end reads ...'
		fastqc --outdir ${path_fastqc_out} ${path_trimm_out}/${filename_trimmed2}
		echo 'Quality checking trimmed data paired end reads completed. Results are stored at' ${path_fastqc_out}
		echo ''
	fi
fi


# Sequence alignment
echo 'Sequence alignment ...'
if [[ ${paired} =~ ^[fF]$ ]] || [[ -z ${filename2} ]]
then
	bwa mem ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed1} > ${path_align_out}/${filename_sam}
elif [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filename2} ]]
then 
	bwa mem ${alignment_settings} ${path_refgenome} ${path_trimm_out}/${filename_trimmed1} ${path_trimm_out}/${filename_trimmed2} > ${path_align_out}/${filename_sam}
fi
echo 'Sequence alignment is completed. Results are stored in' ${path_align_out}/${filename_sam}
echo ''

# Converting sam file to its binary equivalent
echo 'Converting sam to bam ...'
samtools view -b ${path_align_out}/${filename_sam} > ${path_align_out}/${filename_bam}
echo 'Converting sam to bam completed. Results are stored in' ${path_align_out}/${filename_bam}
echo ''

# Quickcheck bam file
echo 'Checking bam file ...'
samtools quickcheck ${path_align_out}/${filename_bam}
echo ''

# Indexing and sorting bam file
if [[ ${sort_and_index} =~ ^[tT]$ ]]
then
	echo 'Indexing bam file ...'
	sambamba-0.7.1-linux-static sort -m 500MB ${path_align_out}/${filename_bam}
	echo 'Indexing completed. Results are stored in' ${path_align_out}
	echo ''
fi


# Transposon mapping
if [[ ${mapping} =~ ^[tT]$ ]]
then
	echo 'Transposon mapping ...'
	cd ~/Documents/Software/python_codes
	python3 ${path_python_codes}transposonmapping_satay.py ${path_align_out}/${filename_sort}
	cd ~/
	echo ''
	echo 'Transposon mapping complete. Results are stored in' ${path_align_out}
	echo ''
fi



if [[ ${delete_sam} =~ ^[tT]$ ]]
then
	echo 'Removing .sam file ...'
	rm ${path_align_out}/${filename_sam}
	echo 'sam file removed.'
fi






### Creating log file
echo ''
echo 'Creating log file ...'
echo ${filename1}	$(date +%F_%T) > ${pathdata}/${filename1%$extension*}'_log.txt'
if [[ ${paired} =~ ^[tT]$ ]] && ! [[ -z ${filename2} ]]
then
	echo 'Paired end reads with paired file:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo ${filename2} >> ${pathdata}/${filename1%$extension*}'_log.txt'
elif [[ ${paired} =~ ^[tT]$ ]] && [[ -z ${filename2} ]]
then
	echo 'Paired end reads with paired reads in same file' >> ${pathdata}/${filename1%$extension*}'_log.txt'
fi

echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Trimming options:' >> ${pathdata}/${filename1%$extension*}'_log.txt'

if [[ ${trimming_software} =~ ^[bB]$ ]]
then
	echo 'BBDuk' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo ${trimming_settings_bbduk} >> ${pathdata}/${filename1%$extension*}'_log.txt'
elif [[ ${trimming_software} =~ ^[tT]$ ]]
then
	echo 'Trimmomatic' >> ${pathdata}/${filename1%$extension*}'_log.txt'
	echo ${trimmomatic_initialization} ${trimming_settings_trimmomatic} >> ${pathdata}/${filename1%$extension*}'_log.txt'
fi

echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Alignment options:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo ${alignment_settings} >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Reference genome used:' ${name_refgenome} >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo '' >> ${pathdata}/${filename1%$extension*}'_log.txt'
echo 'Adapter sequences from adapters.fa:' >> ${pathdata}/${filename1%$extension*}'_log.txt'
cat ${path_bbduk_adapters} >> ${pathdata}/${filename1%$extension*}'_log.txt'





echo 'Processing completed.'



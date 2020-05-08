# LaanLab-SATAY-DataAnalysis
This repository contains codes, data and workflows for data analysis regarding SATAY experiments.

## Projects tab and Issues tab
For comments, ToDo lists and general progress updates on the projects, please notify this in the Projects tab.
Issues and questions about codes can be published in the Issues tab.

## Docs
Files with detailed notes about analysis, datafiles, software and workflows (including an installation guide).

## Data files
Currently the following data files and folder are present:

1. Cerevisiae_EssentialGenes_list_1.txt: This is a list of known essential genes with systematic naming format.
2. Cerevisiae_EssentialGenes_List_2.txt: This is a list of known essential genes with common naming format. Some genes may occur only in one file, so it is recommended to use both files simultaneously to have a complete list of the known essential genes. 
3. Yeast_Protein_Names.txt: This is a list that includes all genes with both naming convention and their ID's. It also include the length of the corresponding proteins in terms of amino acids.
4. S288C_reference_sequence_R64-2-1_20150113.fsa: Reference sequence for wild type cells of *S.Cerevisiae* from the S288C strain.

### test_data
This folder (located in the Data files folder) contains three files that are the result of the processing by the Matlab code of the Kornmann lab (see 'matlab codes for transposon mapping').

1. E-MTAB-4885.WT1.bam.wig: Contains the number of reads per insertion location.
2. E-MTAB-4885.WT1.bam.bed: Contains the number of reads oer insertion location. This file is similar to the .wig file, but in the .bed file the number of reads are given by (reads*20)+100
3. E-MTAB-4885.WT1.bam_pergene.txt: Contain the number of insertions and number of reads per gene. For this an overview is created of most known genes (6603 genes are considered).

## Python scripts and modules
Data analysis created with python are located in the folder for python scripts.
Python functions are placed the python modules folder.

## Matlab codes for transposon mapping
Currently the following Matlab codes and data files are present:

1. tn_and_reads_per_gene.m: This inputs a .bam file and outputs the number of transposons and reads per gene.
2. names.mat: Matlab data file with all genes names. This file is required for running the matlab code 1.
3. yeastGFF.mat: Matlab data file with information about the genes. This file is required for running the matlab code 1.

*Last updated: May 8, 2020*

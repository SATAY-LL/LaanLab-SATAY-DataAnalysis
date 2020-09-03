# LaanLab-SATAY-DataAnalysis

This repository contains codes, data and workflows for data analysis regarding SATAY experiments.

## Projects tab and Issues tab

For comments, ToDo lists and general progress updates on the projects, please notify this in the Projects tab.
Issues and questions about codes can be published in the Issues tab.

## Docs

Files with detailed notes about analysis, datafiles, software and workflows (including an installation guide for the Linux workflow for SATAY analysis).

## Data files

Currently the following data files and folder are present:

1. Cerevisiae_EssentialGenes_list_1.txt: This is a list of known essential genes with systematic naming format.
2. Cerevisiae_EssentialGenes_List_2.txt: This is a list of known essential genes with common naming format. Some genes may occur only in one file, so it is recommended to use both files simultaneously to have a complete list of the known essential genes.
3. Cerevisiae_AllEssentialGenes_list.txt: This is a combination of the previous two files where all the duplicates have been removed and the standard yeast gene naming convention has been applied.
4. Yeast_Protein_Names.txt: This is a list that includes all genes with both naming convention and their ID's. It also include the length of the corresponding proteins in terms of amino acids.
5. S288C_reference_sequence_R64-2-1_20150113.fsa: Reference sequence for wild type cells of *S.Cerevisiae* from the S288C strain.
6. Saccharomyces_cerevisiae.R64-1-1.99.gff3: Information file for all CDS in the *S.Cerevisiae* genome in .gff format.

### test_data

This folder (located in the Data files folder) contains three files that are the result of the processing by the Matlab code of the Kornmann lab (see 'matlab codes for transposon mapping').

1. ... .wig: Contains the number of reads per insertion location in wig format.
2. ... .bed: Contains the number of reads oer insertion location in bed format. This file is similar to the .wig file, but in the .bed file the number of reads are given by (reads*20)+100
3. ... _pergene.txt: Contain the number of insertions and number of reads per gene. For this an overview is created of most known genes (6600 genes are considered).
4. ... _peressential.txt: Contain the number of insertions and number of reads per essential gene. For this an overview is created for the annotated essential genes (1186 genes are considered).
5. ... _pergene_insertions.txt: Contain the exact locations and reads for each insertion in each gene (6600 genes are considered).
6. ... _peressential_insertions.txt: Contain the exact locations and reads for each insertion in each annotated essential gene (6600 genes are considered).

## Python scripts and python modules

Data analysis created with python are located in the folder for python scripts.
Python functions are placed the python modules folder.
Each python code should include information on how and when to use it.
Some codes also come a python notebook format using the .ipynb extension.

## Matlab codes for transposon mapping

Currently the following Matlab codes and data files are present:

1. tn_and_reads_per_gene.m: This inputs a .bam file and outputs the number of transposons and reads per gene.
2. names.mat: Matlab data file with all genes names. This file is required for running the matlab code 1.
3. yeastGFF.mat: Matlab data file with information about the genes. This file is required for running the matlab code 1.

## Python codes for transposon mapping

This python is used in the Linux workflow for analyzing SATAY data.
The file transposonmapping_satay.py includes the actual code and the folder python_modules include some modules necessary to run it.
This folder also include processing_workflow shell script that is used to call all software tools (including this python script) within Linux.
For setting up a virtual environment including all software tools required for running this workflow, see the installation guide in the docs folder on this Github page.
Note: The python code only works in Linux operating systems.

*Last updated: August 3, 2020*

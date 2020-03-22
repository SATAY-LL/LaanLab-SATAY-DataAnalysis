# LaanLab-SATAY-DataAnalysis
This repository contains codes, data and workflows for data analysis regarding SATAY experiments.

## Projects tab and Issues tab
For comments, ToDo lists and general progress updates on the projects, please notify this in the Projects tab.
Issues about codes can be published in the Issues tab.

## Data files
Currently the following data files are present:

1. Cerevisiae_EssentialGenes_list_1.txt: This is a list of known essential genes with systematic naming format.
2. Cerevisiae_EssentialGenes_List_2.txt: This is a list of known essential genes with common naming format. Some genes may occur only in one file, so it is recommended to use both files simultaneously to have a complete list of the known essential genes. 
3. Yeast_Protein_Names.txt: This is a list that includes all genes with both naming convention and their ID's. It also include the length of the corresponding proteins in terms of amino acids.
4. S288C_reference_sequence_R64-2-1_20150113.fsa: Reference sequence for wild type cells of *S.Cerevisiae* from the S288C strain.
5. E-MTAB-4885.WT1.bam_pergene.txt: This file is an example of how the result looks like after processing and can be used as a test file. It contains a list of all analyzed genes with the respective read and transposon counts. This file can be input in the statistics_pergene.py code.

## Python notebooks
Currently the following notebooks are present:

1. statistics_pergene.ipynb & statistics_pergene.py: This file reads all genes from the text file that is given as output from the Matlab code of Benoit. It takes the number of reads and transposons per gene normalized for the length of the gene and determine statistical values for this.

## Matlab codes
Currently the following Matlab codes and data files are present:

1. tn_and_reads_per_gene.m: This inputs a .bam file and outputs the number of transposons and reads per gene. The output of this code can be used for python notebook 1.
2. names.mat: Matlab data file with all genes names. This file is required for running the matlab code 1.
3. yeastGFF.mat: Matlab data file with information about the genes. This file is required for running the matlab code 1.

*Last updated: March 22, 2020*

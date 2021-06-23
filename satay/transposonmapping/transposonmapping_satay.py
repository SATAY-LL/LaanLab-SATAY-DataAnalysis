#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a tool developed for analysing transposon insertions for experiments using SAturated Transposon Analysis in Yeast (SATAY).
This python code contains one function called transposonmapper().
For more information about this code and the project, see github.com/Gregory94/LaanLab-SATAY-DataAnalysis

This code is based on the Matlab code created by the Kornmann lab which is available at: sites.google.com/site/satayusers/

__Author__ = Gregory van Beek. LaanLab, department of Bionanoscience, Delft University of Technology
__version__ = 1.5
__Date last update__ = 2021-01-11

Version history:
    1.1; Added code for creating two text files for storing insertion locations per gene and per essential gene [2020-07-27]
    1.2; Improved searching algorithm for essential genes [2020-08-06]
    1.3; Load file containing all essential genes so that a search for essential genes in multiple file is not needed anymore. This file is created using Create_EssentialGenes_list.py located in the same directory as this code [2020-08-07]
    1.4; Fixed bug where the gene position and transposon insertion location did not start at zero for each chromosome, causing confusing values to be stored in the _pergene_insertions.txt and _peressential_insertions.txt files [2020-08-09]
    1.5; Added functionality to handle all possible sam flags in the alignment file (bam-file) instead of only flag=0 or flag=16. This is needed for the function to handle paired-end sequencing data [2021-01-11]
"""

import pysam

# Local imports
from .python_modules import chromosomename_roman_to_arabic
from .python_modules import get_chromosome_names
from .python_modules import get_sequence_length
from .python_modules import get_reads
from .python_modules import add_chromosome_length
from .python_modules import add_chromosome_length_inserts
from .python_modules import get_insertions_and_reads

from .exporting import (
    save_as_bed,
    save_per_gene,
    save_per_gene_insertions,
    save_per_essential_insertions,
    export_as_wig,
)

from .importing import load_default_files, read_genes


def transposonmapper(bamfile, gff_file=None, essential_file=None, gene_name_file=None):
    """
    This function is created for analysis of SATAY data using the species Saccharomyces Cerevisiae.
    It outputs the following files that store information regarding the location of all insertions:
        - .bed-file: Includes all individual basepair locations of the whole genome where at least one transposon has been mapped and the number of insertions for each locations (the number of reads) according to the Browser Extensible Data (bed) format.
                    A distinction is made between reads that had a different reading orientation during sequencing. The number of reads are stored using the equation #reads*20+100 (e.g. 2 reads is stored as 140).
        - .wig-file: Includes all individual basepair locations of the whole genome where at least one transposon has been mapped and the number of insertions for each locations (the number of reads) according to the Wiggle (wig) format.
                    In this file no distinction is made between reads that had a different reading orientation during sequencing. The number of reads are stored as the absolute count.
        - _pergene.txt-file: Includes all genes (currently 6600) with the total number of insertions and number of reads within the genomic region of the gene.
        - _peressential.txt-file: Includes all annotated essential genes (currently 1186) with the total number of insertions and number of reads within the genomic region of the gene.
        - _pergene_insertions.txt-file: Includes all genes with their genomic location (i.e. chromosome number, start and end position) and the locations of all insertions within the gene location. It also include the number number of reads per insertions.
        - _peressential_insertions.txt-file: Includes all essential genes with their genomic location (i.e. chromosome number, start and end position) and the locations of all insertions within the gene location. It also include the number number of reads per insertions.
          (note that in the latter two files, the genomic locations are continous, for example chromosome II does not start at 0, but at 'length chromosome I + 1' etc.).
    The output files are saved at the location of the input file using the same name as the input file, but with the corresponding extension.
    
    The function assumes that the reads are already aligned to a reference genome.
    The input data should be a .bam-file and the location where the .bam-file is stored should also contain an index file (.bam.bai-file, which for example can be created using sambamba).
    This function takes the following inputs:
        - bamfile [required]: Path to the bamfile. This location should also contain the .bam.bai index file (does not need to be input in this function).
        - gfffile [optional]: Path to a .gff-file including all gene information (e.g. downloaded from SGD). Default file is 'Saccharomyces_cerevisiae.R64-1-1.99.gff3'.
        - essentialfiles [optional]: Path to a .txt file containing a list all essential genes. Every line should consist of a single essential gene and the file should have one header line. Ideally this file is created using 'Create_EssentialGenes_list.py'. Default file is 'Cerevisiae_AllEssentialGenes_List.txt'.
        - genenamesfile [optional]: Path to text file that includes aliases for all genes. Default file is 'Yeast_Protein_Names.txt'.
    When the arguments for the optional files are not given, the files are used that are stored at the following location:
        "path_current_pythonscript/../data_files"
    The function uses the pysam package for handling bam files (see pysam.readthedocs.io/en/latest/index.html) and therefore this function only runs on Linux systems with SAMTools installed.
    """

    # If necessary, load default files
    gff_file, essential_file, gene_name_file = load_default_files(
        gff_file, essential_file, gene_name_file
    )

    # Verify presence of files
    data_files = {
        "bam": bamfile,
        "gff3": gff_file,
        "essentials": essential_file,
        "gene_names": gene_name_file,
    }

    for filetype, file_path in data_files.items():
        assert file_path, f"{filetype} not found at {file_path}"

    # Read bam file
    bam = pysam.AlignmentFile(bamfile, "rb")

    # Get names of all chromosomes as stored in the bam file
    ref_tid = get_chromosome_names(bam)
    ref_names = list(ref_tid.keys())

    # Convert chromosome names in data file to roman numerals
    ref_romannums = chromosomename_roman_to_arabic()[1]
    ref_tid_roman = {key: value for key, value in zip(ref_romannums, ref_tid)}

    # Get sequence lengths of all chromosomes
    chr_lengths, chr_lengths_cumsum = get_sequence_length(bam)

    # Get all reads within a specified genomic region
    readnumb_array, tncoordinates_array, tncoordinatescopy_array = get_reads(bam)

    # Read files for all genes and all essential genes
    print("Getting coordinates of all genes ...")
    gene_coordinates, essential_coordinates, aliases_designation = read_genes(
        gff_file, essential_file, gene_name_file
    )

    #%% CONCATENATE ALL CHROMOSOMES

    # For each insertion location, add the length of all previous chromosomes
    tncoordinatescopy_array = add_chromosome_length_inserts(
        tncoordinatescopy_array, ref_names, chr_lengths
    )

    # For each gene location, add the length of all previous chromosomes
    gene_coordinates = add_chromosome_length(
        gene_coordinates, chr_lengths_cumsum, ref_tid_roman
    )

    # For each essential gene location, add the length of all previous chromosomes
    essential_coordinates = add_chromosome_length(
        essential_coordinates, chr_lengths_cumsum, ref_tid_roman
    )

    # GET NUMBER OF TRANSPOSONS AND READS PER GENE
    print("Get number of insertions and reads per gene ...")

    # All genes
    tn_per_gene, reads_per_gene, tn_coordinates_per_gene = get_insertions_and_reads(
        gene_coordinates, tncoordinatescopy_array, readnumb_array
    )

    # Only essential genes
    (
        tn_per_essential,
        reads_per_essential,
        tn_coordinates_per_essential,
    ) = get_insertions_and_reads(
        essential_coordinates, tncoordinatescopy_array, readnumb_array
    )

    # CREATE BED FILE
    bedfile = bamfile + ".bed"
    print("Writing bed file at: ", bedfile)
    print("")
    save_as_bed(bedfile, tncoordinates_array, ref_tid, readnumb_array)

    # CREATE TEXT FILE WITH TRANSPOSONS AND READS PER GENE
    # NOTE THAT THE TRANSPOSON WITH THE HIGHEST READ COUNT IS IGNORED.
    # E.G. IF THIS FILE IS COMPARED WITH THE _PERGENE_INSERTIONS.TXT FILE THE READS DON'T ADD UP (SEE https://groups.google.com/forum/#!category-topic/satayusers/bioinformatics/uaTpKsmgU6Q)
    # TOO REMOVE THIS HACK, CHANGE THE INITIALIZATION OF THE VARIABLE readpergene
    per_gene_file = bamfile + "_pergene.txt"
    print("Writing pergene.txt file at: ", per_gene_file)
    print("")

    save_per_gene(per_gene_file, tn_per_gene, reads_per_gene, aliases_designation)

    # CREATE TEXT FILE TRANSPOSONS AND READS PER ESSENTIAL GENE
    per_essential_file = bamfile + "_peressential.txt"
    print("Writing peressential.txt file at: ", per_essential_file)
    print("")
    save_per_gene(
        per_essential_file, tn_per_essential, reads_per_essential, aliases_designation
    )

    # CREATE TEXT FILE WITH LOCATION OF INSERTIONS AND READS PER GENE
    per_gene_insertions_file = bamfile + "_pergene_insertions.txt"
    print("Witing pergene_insertions.txt file at: ", per_gene_insertions_file)
    print("")

    save_per_gene_insertions(
        per_gene_insertions_file,
        tn_coordinates_per_gene,
        gene_coordinates,
        chr_lengths_cumsum,
        ref_tid_roman,
        aliases_designation,
    )

    # CREATE TEXT FILE WITH LOCATION OF INSERTIONS AND READS PER ESSENTIAL GENE
    per_essential_insertions_file = bamfile + "_peressential_insertions.txt"
    print(
        "Writing peressential_insertions.txt file at: ", per_essential_insertions_file
    )
    print("")

    save_per_essential_insertions(
        per_essential_insertions_file,
        tn_coordinates_per_essential,
        gene_coordinates,
        chr_lengths_cumsum,
        ref_tid_roman,
        aliases_designation,
    )

    # ADD INSERTIONS AT SAME LOCATION BUT WITH DIFFERENT ORIENTATIONS TOGETHER (FOR STORING IN WIG-FILE)
    wigfile = bamfile + ".wig"
    print("Writing wig file at: ", wigfile)
    print("")

    export_as_wig(wigfile, ref_names, readnumb_array, tncoordinates_array, ref_tid)


#%%
if __name__ == "__main__":
    bamfile = "satay/data_files/files4test/SRR062634.filt_trimmed.sorted.bam"
    transposonmapper(bamfile=bamfile)


# SATAY analysis workflow

[![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)

## Overview

This workflow is created for processing sequencing data for SAturated Transposon Analysis in Yeast for Saccharomyces Cerevisiae.
For more information about this project, see our [JupyterBook](https://leilaicruz.github.io/SATAY-jupyter-book/Introduction.html).
For more information regarding SATAY, see [the satay user website](https://sites.google.com/site/satayusers/) created by the Kornmann-lab.

## Features

Requires input sequencing data in fastq format.
It can perform the following tasks:
- sequence trimming
- quality checking raw and trimmed fastq files
- sequence alignment with reference genome (S288C Cerevisiae genome)
- quality checking, indexing and sorting of alignment
- transposon mapping

Output files indicate the location of transposon insertions and the number of reads at those locations.
This is presented in both .bed and .wig format.
Also a list of genes is generated where the number or distribution of insertions and reads is presented per (essential) gene.

## Running the project

...

## Dependencies

- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.9 or later
- [BBMap](https://sourceforge.net/projects/bbmap/) v38.87 or later
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.39 or later
- [BWA](https://sourceforge.net/projects/bio-bwa/) v0.7.17 or later
- [SAMTools](http://www.htslib.org/download/) v1.10 or later
- [BCFTools](http://www.htslib.org/download/) v1.10.2-3 or later
- [Sambamba](https://github.com/biod/sambamba/releases) v0.7.1 or later
- Python v3.7 or later
  - [pysam](https://anaconda.org/bioconda/pysam) v0.16.0.1 or later
  - [numpy](https://anaconda.org/anaconda/numpy) v1.19.2 or later

## ToDo

- [x] Create GUI
- [] Test working environment in Linux

## Contributors
[LaanLab. Department of BioNanoScience, Delft University of Technology](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/bionanoscience/research/research-labs/liedewij-laan-lab/research-projects/evolvability-and-modularity-of-essential-functions-in-budding-yeast)
- Leila Inigo de la Cruz
- Enzo Kingma
- Gregory van Beek

*Last updated: February 15, 2021*

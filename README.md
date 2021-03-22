# SATAY analysis workflow

## Overview

This workflow is created for processing sequencing data for SAturated Transposon Analysis in Yeast (SATAY) for Saccharomyces Cerevisiae.
It performs the steps from raw sequencing data until the transposon mapping that outputs files containing all insertion sites combined with the number of reads.
For more information about this project, see our [JupyterBook](https://leilaicruz.github.io/SATAY-jupyter-book/Introduction.html).
For more information regarding SATAY, see [the satay user website](https://sites.google.com/site/satayusers/) created by the Kornmann-lab.

For a complete guide how to use it, see the [documentation](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/documentation/documentation_satay.md).

## Features

The workflow requires input sequencing data in fastq format.
It can perform the following tasks:

- sequence trimming
- quality checking raw and trimmed fastq files
- sequence alignment with reference genome (S288C Cerevisiae genome)
- quality checking bam files, indexing and sorting
- transposon mapping

The output files indicate the location of transposon insertions and the number of reads at those locations.
This is presented in both .bed and .wig format.
Also a list of genes is generated where the number and distribution of insertions and reads is presented per (essential) gene.

## Running the project

The main workflow is called [satay.sh](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/satay.sh) and runs in Linux.
Before starting the workflow on your Linux machine, you might need to change some paths in the script.
Start the workflow using the command `bash satay.sh`.

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

## Contributors

[LaanLab. Department of BioNanoScience, Delft University of Technology](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/bionanoscience/research/research-labs/liedewij-laan-lab/research-projects/evolvability-and-modularity-of-essential-functions-in-budding-yeast)

- Leila Inigo de la Cruz
- Enzo Kingma
- Gregory van Beek

## License

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
This work is licensed under Apache 2.0 .
The 2.0 version of the Apache License, approved by the ASF in 2004, helps us achieve our goal of providing reliable and long-lived software products through collaborative open source software development.

*Last updated: March 3, 2021*
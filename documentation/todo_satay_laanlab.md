# SATAY processing - To Do

> Date last update: 26-03-2021
>
> Author: Gregory van Beek

- [ ] Make a plan how to define trimming and alignment settings.
Currently the trimming and processing options are chosen by trial-and-error.
Each dataset requires unique processing and therefore no fixed settings for the trimming and alignment can be determined, but potentially some guidelines can be created that helps for the processing.

- [ ] The code is relying on software third party software tools (e.g. for the trimming, alignment and quality checking).
This currently requires to be installed manually and in the [workflow](https://github.com/leilaicruz/LaanLab-SATAY-DataAnalysis/blob/master/satay.sh), paths need to be set to some these tools (in satay.sh, the paths are all defined in the beginning of the script in the `DEFINE PATHS` section). This is not convenient for other users and is prone to errors.
Therefore a package need to be created that installs the third party software tools at a fixed location relative to the workflow.
Then the paths doesn't have to set manually by the user and makes it easier to set up and use.

- [ ] Implementing unit testing in the Github repository. This allows for automatic testing of the codes after changes have been made in the repository and should ensure the outcome is still correct.

- [ ] For the transposon mapping in the pipeline [a custom python script](https://github.com/leilaicruz/LaanLab-SATAY-DataAnalysis/blob/master/python_transposonmapping/transposonmapping_satay.py) is created that inputs a bam file and outputs lists of insertion locations and the corresponding number of reads. This script is based on the [matlab code](https://sites.google.com/site/satayusers/complete-protocol/bioinformatics-analysis/matlab-script) from the Kornmann lab.
The matlab code contained some bugs that were solved in the python script (see the [satay forum](https://groups.google.com/g/satayusers/search?q=matlab)). Therefore some small differences may be present between the matlab processing and the python processing. But recently, some more differences were found that could not be explained by the matlab bugs.
The exact cause is yet unclear, but has maybe to do with reading the samflag in the python script.

- [ ] Currently [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) is used for the alignment of the reads.
This tools work fine, but the Kornmann lab (and many other labs) use [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml).
This is already installed on the Linux machine, but needs to be implemented in the [workflow](https://github.com/leilaicruz/LaanLab-SATAY-DataAnalysis/blob/master/satay.sh).
The current aligner is implemented in the section `SEQUENCE ALIGNMENT`.
Note to set the option for paired-end and single-end alignment for bowtie (defined in the variabele `${paired}`).

- [ ] Automated documentation based on the help texts from the python scripts.

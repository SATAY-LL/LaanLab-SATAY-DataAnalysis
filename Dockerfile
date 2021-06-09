FROM continuumio/miniconda3

SHELL [ "/bin/bash", "--login", "-c" ]

RUN apt-get update -q && apt-get install yad xdg-utils -q -y && apt-get clean 

RUN conda install --quiet --yes --freeze-installed \ 
    -c conda-forge -c bioconda -c agbiome \
    fastqc \
    trimmomatic \
    bbtools \
    bwa \
    samtools \
    bcftools \
    sambamba \
    pysam \
    numpy \
    && conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete 

ENV adapters=/opt/conda/bbtools/lib/resources/adapters.fa
ENV bbduk=/opt/conda/bbtools/lib/bbduk.sh

# COPY . /opt/satay

RUN git clone https://github.com/SATAY-LL/LaanLab-SATAY-DataAnalysis.git /opt/satay

CMD bash /opt/satay/satay.sh

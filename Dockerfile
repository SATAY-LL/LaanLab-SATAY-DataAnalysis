FROM continuumio/miniconda3

SHELL [ "/bin/bash", "--login", "-c" ]

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

# COPY . /opt/satay

RUN git clone https://github.com/SATAY-LL/LaanLab-SATAY-DataAnalysis.git /opt/satay

CMD bash /opt/satay/satay.sh -h

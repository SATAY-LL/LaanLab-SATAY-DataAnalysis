from satay.transposonmapping.properties.get_gene_position import gene_position
from satay.transposonmapping.properties.gene_aliases import gene_aliases


def read_genes(gff_file, essentials_file, gene_names_file):
    """
    This function reads the useful information inside the gff_file, essentials_file and
    gene_names_file. 
    For the gff_file and essentials_file extracts the gene coordinates , specifying the chromosome, start ,end 
    and direction. 
    For the gene_names_files it translates the systematic name into the standard name.
    
    Usage: 

     [gene_coordinates, essential_coordinates, aliases_designation] = 
            read_genes(gff_file, essentials_file, gene_names_file)
    
    Inputs:
        - gff_file : Annotated genome from Saccharomyces cerevisiae (baker's yeast) (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz)
        - essentials_file:  Essentials genes annotated from yeast , written using the systematic name standard , all in one column
        - gene_names_file:  This documents lists all the Saccharomyces cerevisiae S288c entries present in
this release of UniProtKB/Swiss-Prot. 
            Yeast (Saccharomyces cerevisiae): entries, gene names and cross-references to SGD. 
             Release:     2021_01 of 10-Feb-2021
            
    Outputs:
        gene_coordinates : a dict specifying for each gene the chromosome number the gene
        belongs to, the start gene coordinate, the end gene coordinate and the strand direction ('+' or '-'). 
        essential_coordinates: a dict specifying for each annotated essential gene the chromosome number the gene
        belongs to, the start gene coordinate, the end gene coordinate and the strand direction ('+' or '-'). 
        aliases_designation: a dict that for each systematic gene name specify the standard gene name.  

    """

    # Get gene position
    gene_coordinates = gene_position(gff_file)  #'YAL069W' | ['I', 335, 649], ...

    # Get all annotated essential genes
    essential_coordinates = {}
    with open(essentials_file, "r") as f:
        genes = f.readlines()[1:]
        for gene in genes:
            name = gene.strip("\n")
            essential_coordinates[name] = gene_coordinates.get(name).copy()

    # Get aliases of all genes
    aliases_designation = gene_aliases(gene_names_file)[0]  #'YMR056C' \ ['AAC1'], ...

    return gene_coordinates, essential_coordinates, aliases_designation

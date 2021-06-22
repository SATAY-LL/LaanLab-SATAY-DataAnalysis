from .chromosome_and_gene_positions import gene_position
from .gene_names import gene_aliases

def read_genes(gff_file, essentials_file, gene_names_file):
    """

    Syntax: [gene_coordinates, essential_coordinates, aliases_designation] = 
            read_genes(gff_file, essentials_file, gene_names_file)

    """
    
    # Get gene position
    gene_coordinates = gene_position(gff_file) #'YAL069W' | ['I', 335, 649], ...

    # Get all annotated essential genes
    essential_coordinates = {}
    with open(essentials_file, 'r') as f:
        genes = f.readlines()[1:]
        for gene in genes:
            name = gene.strip('\n')
            essential_coordinates[name] = gene_coordinates.get(name).copy()


    # Get aliases of all genes
    aliases_designation = gene_aliases(gene_names_file)[0] #'YMR056C' \ ['AAC1'], ...

    return gene_coordinates, essential_coordinates, aliases_designation

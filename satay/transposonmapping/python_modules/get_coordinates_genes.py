import os

from satay.transposonmapping.python_modules.chromosome_and_gene_positions import (
    gene_position,
)
from satay.transposonmapping.python_modules.gene_names import gene_aliases


def get_coordinates_genes(path, data_files):
    """
    Get coordinates of all genes

    Parameters
    ----------
    path : str
    data_files : dict

    Returns
    -------
    essential_coordinates : dict    
    aliases_designation : dict

    """

    essential_coordinates = {}

    # Get position genes
    if "gff3" in data_files:
        file_path = os.path.join(path, data_files["gff3"])
        gene_coordinates = gene_position(file_path)
    else:
        raise ValueError("gff3 type not found in data")

    # Get all annotated essential genes
    if "essential_genes" in data_files:
        file_path = os.path.join(path, data_files["essential_genes"])
        with open(file_path, "r") as f:
            genes = f.readlines()[1:]
            for gene in genes:
                name = gene.strip("\n")
                essential_coordinates[name] = gene_coordinates.get(name).copy()
    else:
        raise ValueError("essential_genes not found in data")

    # Get aliases of all genes
    if "gene_names" in data_files:
        file_path = os.path.join(path, "Yeast_Protein_Names.txt")
        aliases_designation = gene_aliases(file_path)[0]  #'YMR056C' \ ['AAC1'], ...
    else:
        raise ValueError("gene_names not found in data")

    return essential_coordinates, aliases_designation

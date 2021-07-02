import os
import pkg_resources


def load_default_files(gff_file, essentials_file, gene_names_file):
    """

    """

    default_path = pkg_resources.resource_filename("satay", "data_files/")
    if gff_file is None:
        gff_file = os.path.join(
            default_path, "Saccharomyces_cerevisiae.R64-1-1.99.gff3"
        )
    if essentials_file is None:
        essentials_file = os.path.join(
            default_path, "Cerevisiae_AllEssentialGenes_List.txt"
        )
    if gene_names_file is None:
        gene_names_file = os.path.join(default_path, "Yeast_Protein_Names.txt")


    return gff_file, essentials_file, gene_names_file


def load_sgd_tab(sgd_features_file):
    """
    
    """
    default_path = pkg_resources.resource_filename("satay", "data_files/")
    if sgd_features_file is None: 

        sgd_features_file = os.path.join(default_path,'SGD_features.tab')

    return sgd_features_file


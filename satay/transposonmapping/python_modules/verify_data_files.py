import os


def verify_data_files(path="", data_files=None):
    """
    Verify presence of the following essential data files
    - bam
    - gff3
    - essential genes
    - gene names

    Parameters
    ----------
    path : str
        Folder containing all data files.
    data_files : dict
        Dictionary containing filetype and filename as key:value

    TODO: put data_files in a separate config file

    """

    # Verify existence data folder
    assert os.path.isdir(path), f"{path} was not found."

    if data_files is None:
        data_files = {
            # "bamfile": "WT_merged-techrep-a_techrep-b_trimmed.sorted.bam",
            "gfffile": "Saccharomyces_cerevisiae.R64-1-1.99.gff3",
            "essential_genes": "Cerevisiae_AllEssentialGenes_List.txt",
            "gene_names": "Yeast_Protein_Names.txt",
        }

    # Check presence data files
    for filetype, file_name in data_files.items():
        file_path = os.path.join(path, file_name)
        assert os.path.isfile(file_path), f"{filetype} not found at {file_path}"

import os
import pkg_resources

class Files:
    """ Container for transposon mapping datafiles


    """

    def __init__(
        self, bam_file=None, gff_file=None, essentials_file=None, gene_names_file=None,
    ):
        self.default_path = pkg_resources.resource_filename("satay", "data_files/")
        
        self.bam_file = bam_file             
        self.gff_file = gff_file 
        self.essentials_file = essentials_file
        self.gene_names_file = gene_names_file
        
        self._set_default_files()
        self.verify_data_files()

    def _set_default_files(self):        
        if self.gff_file is None:
            self.gff_file = os.path.join(self.default_path, "Saccharomyces_cerevisiae.R64-1-1.99.gff3")
        if self.essentials_file is None:
            self.essentials_file = os.path.join(self.default_path, "Cerevisiae_AllEssentialGenes_List.txt")
        if self.gene_names_file is None:
            self.gene_names_file = os.path.join(self.default_path, "Yeast_Protein_Names.txt")

    def verify_data_files(self):        
        data_files = {
            "bam": self.bam_file,
            "gff3": self.gff_file,
            "essentials": self.essentials_file,
            "gene_names": self.gene_names_file,
        }
        for filetype, file_path in data_files.items():
            assert file_path, f"{filetype} not found at {file_path}"

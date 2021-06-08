import pkg_resources

# import pysam

from .python_modules import verify_data_files


class Data:
    """ Container for transposon mapping datafiles


    """

    def __init__(
        self, bam_file="", gff_file="", essentials_file="", gene_names_file="",
    ):
        self.bam_file = (
            bam_file if bam_file else "WT_merged-techrep-a_techrep-b_trimmed.sorted.bam"
        )
        self.gff_file = (
            gff_file if gff_file else "Saccharomyces_cerevisiae.R64-1-1.99.gff3"
        )
        self.essentials_file = (
            essentials_file
            if essentials_file
            else "Cerevisiae_AllEssentialGenes_List.txt"
        )
        self.gene_names_file = (
            gene_names_file if gene_names_file else "Yeast_Protein_Names.txt"
        )
        self.path = pkg_resources.resource_filename("satay", "data_files/")

        self.verify_data_files()
        self.get_bam()

    def verify_data_files(self):
        data_files = {
            # "bam": bam_file,
            "gff3": self.gff_file,
            "essentials": self.essentials_file,
            "gene_names": self.gene_names_file,
        }
        verify_data_files(path=self.path, data_files=data_files)

    def get_bam(self):
        # self.bam = pysam.AlignmentFile(self.bamfile, "rb")
        pass

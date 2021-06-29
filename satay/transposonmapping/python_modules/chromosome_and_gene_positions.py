from ..utils import chromosomename_roman_to_arabic


"""This module includes three functions, 'chromosome_position', 'chromosomename_roman_to_arabic' and 'gene_position'.
Except for 'chromosomename_roman_to_arabic', the functions require an input with the full path to the file 'Saccharomyces_cerevisiae.R64-1-1.99.gff3' (which can be downloaded from https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index).
If no input is given, the program automatically searches for the file at in the current directory.
The function 'chromosomename_roman_to_arabic' does not take any input. It returns two conversion lists for numbers 1 to 16 to get from arabic to roman numerals or vice versa, respectively.
"""

#%%
def chromosome_position(gff_file):
    """Get the start and end position of each chromosome and determine their respective length.
    Input is a .gff file downloaded from https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index
    Output are three dictionaries for length, start and end position. All dictionaries have keys representing the chromosome number in roman numerals.
    To get all dictionaries, use: 'a,b,c, chromosome_and_gene_position.chromosome_position()'.
    'a' = chromosome length
    'b' = chromosome start position
    'c' = chromosome end position
    """

    # Get roman to arabic dictionary
    roman_to_arabic = chromosomename_roman_to_arabic()[1]

    # GET END POSITIONS OF THE CHROMOSOMES FROM THE GFF FILE AND STORE THEM IN A DICTIONARY
    chr_length_dict = {}
    chr_length_list = []
    with open(gff_file) as f:
        line_counter = 0
        next(f)
        while line_counter < 17:
            lines = f.readline()
            chr_line_list = lines.strip("\n").replace(" ", "\t").split("\t")
            chr_number = chr_line_list[3]
            # if chr_number != 'Mito':
            chr_length = int(chr_line_list[5])
            chr_length_list.append(chr_length)
            chr_length_dict[chr_number] = chr_length
            line_counter += 1

    # DETERMINE START AND END POSITION OF EACH OF THE CHROMOSOMES
    chr_start_pos_dict = {}
    chr_end_pos_dict = {}
    counter = 0
    for roman in roman_to_arabic.keys():
        chr_start_pos_dict[roman] = sum(chr_length_list[:counter]) + 1
        chr_end_pos_dict[roman] = sum(chr_length_list[: counter + 1])
        if roman == "I":
            chr_start_pos_dict[roman] = 1
        counter += 1

    return (chr_length_dict, chr_start_pos_dict, chr_end_pos_dict)


#%%
def gene_position(gff_file, get_dict=True):
    """Get the start and end position of each gene and determine their respective length.
    Input is a .gff file downloaded from https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index
    Output is a dictionary that includes all gene names as keys. The values are lists with four inputs.
    The first is the chromosome number the gene belong to, the second is the start position, the third is the end position of the gene in terms of basepairs, the fourth is the reading orientation of the gene.
    The reading orientation is indicated with a '+' (forward reading) or '-' (reverse reading).
    The get_dict by defulat sets that the output should be given as a dictionary with keys the different genes and the values a list of the different parameters.
    When the get_dict is set to False, the code returns all the values as individual lists.
    """

    if get_dict == True:
        gene_pos_dict = {}
        with open(gff_file) as f:
            for line in f:
                line_list = line.split("\t")
                if len(line_list) > 2:
                    if line_list[2] == "gene":
                        gene_chr = line_list[0]
                        gene_start = line_list[3]
                        gene_end = line_list[4]
                        gene_orien = line_list[6]
                        gene_position = [
                            gene_chr,
                            int(gene_start),
                            int(gene_end),
                            gene_orien,
                        ]

                        gene_name_string = line_list[8].split(";")[0]
                        gene_name = gene_name_string.split(":")[1]

                        gene_pos_dict[gene_name] = gene_position

        return gene_pos_dict

    else:
        gene_chr = []
        gene_start = []
        gene_end = []
        gene_orien = []
        gene_name = []
        with open(gff_file) as f:
            for line in f:
                line_list = line.split("\t")
                if len(line_list) > 2:
                    if line_list[2] == "gene":
                        gene_chr.append(line_list[0])
                        gene_start.append(int(line_list[3]))
                        gene_end.append(int(line_list[4]))
                        gene_orien.append(line_list[6])

                        gene_name_string = line_list[8].split(";")[0]
                        gene_name.append(gene_name_string.split(":")[1])

        return (gene_name, gene_chr, gene_start, gene_end, gene_orien)


#%%
if __name__ == "__main__":
    chromosome_position()

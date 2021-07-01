from ..utils import chromosomename_roman_to_arabic

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

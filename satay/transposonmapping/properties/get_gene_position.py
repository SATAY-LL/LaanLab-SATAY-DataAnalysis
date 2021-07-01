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

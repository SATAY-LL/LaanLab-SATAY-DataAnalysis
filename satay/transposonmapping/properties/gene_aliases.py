def gene_aliases(gene_information_file=None):
    """Create three dictionaries containing aliases for genes
    Input is the path to 'Protein_Names.txt' file downloaded from https://www.uniprot.org/docs/yeast.
    If no input is given the file is automatically searched for at thisscriptlocation/../Data_Files/Yeast_Protein_Names.txt.
    Output is three dictionaries:
    aliases_designation_dict = gene aliases for common names (e.g. Bem1 and Sro1)
    aliases_sgd_dict = gene aliases for the search names in SGD (e.g. Bem1 and S000000404)
    aliases_swissprot_dict = gene aliases for Swiss Prot.
    The keys of the dictionaries are the systematic names of the genes (e.g. YBR200W for Bem1)
    
    Search through lists to get corresponding key:
        [key for key, val in aliases.items() if 'TFC3' in val]
    """

    if gene_information_file == None:
        import os

        file_dirname = os.path.dirname(os.path.abspath("__file__"))
        gene_information_file = os.path.join(
            file_dirname, "..", "..", "data_files", "Yeast_Protein_Names.txt"
        )

    aliases_designation_dict = {}
    aliases_sgd_dict = {}
    aliases_swissprot_dict = {}

    with open(gene_information_file) as f:
        lines = f.readlines()

        for i in range(
            58, len(lines) - 6
        ):  # THE GENES START AT LINE 58 AND STOP 6 LINES BEFORE THE END OF THE FILE.
            l = lines[i]
            l_list = l.split()

            alias_counter = 0
            gene_designation_list = []
            for names in l_list:
                if names.endswith(";"):
                    gene_designation_list.append(names.strip(";"))
                    alias_counter += 1
            gene_designation_list.append(l_list[alias_counter])

            oln_name = l_list[alias_counter + 1]
            sgd_name = l_list[alias_counter + 4]
            swissprot_name = l_list[alias_counter + 2]
            if oln_name == "GAG" or oln_name == "POL":
                oln_name = l_list[alias_counter + 2]
                sgd_name = l_list[alias_counter + 5]
                swissprot_name = l_list[alias_counter + 3]

            aliases_designation_dict[oln_name] = gene_designation_list
            aliases_sgd_dict[oln_name] = sgd_name
            aliases_swissprot_dict[oln_name] = swissprot_name

    return (aliases_designation_dict, aliases_sgd_dict, aliases_swissprot_dict)


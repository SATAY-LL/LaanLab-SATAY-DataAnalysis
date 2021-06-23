def save_per_gene(filename, tn_per_gene, reads_per_gene, aliases_designation):
    """
    Create text file with transposons and reads per gene
 
    NOTE THAT THE TRANSPOSON WITH THE HIGHEST READ COUNT IS IGNORED.
    E.G. IF THIS FILE IS COMPARED WITH THE _PERGENE_INSERTIONS.TXT FILE THE 
    READS DON'T ADD UP (SEE https://groups.google.com/forum/#!category-topic/satayusers/bioinformatics/uaTpKsmgU6Q)
    TOO REMOVE THIS HACK, CHANGE THE INITIALIZATION OF THE VARIABLE readpergene
 
    """

    with open(filename, "w") as f:

        f.write("Gene name\tNumber of transposons per gene\tNumber of reads per gene\n")

        for gene in tn_per_gene:
            tnpergene = tn_per_gene[gene]
            readpergene = reads_per_gene[gene]
            if gene in aliases_designation:
                gene_alias = aliases_designation.get(gene)[0]
            else:
                gene_alias = gene
            f.write(gene_alias + "\t" + str(tnpergene) + "\t" + str(readpergene) + "\n")


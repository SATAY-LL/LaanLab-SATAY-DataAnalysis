def save_per_gene_insertions(
    filename,
    tn_coordinates,
    gene_coordinates,
    chr_lengths_cumsum,
    ref_tid_roman,
    aliases_designation,
):
    """

    
    """

    with open(filename, "w") as f:

        f.write(
            "Gene name\tChromosome\tStart location\tEnd location\tInsertion locations\tReads per insertion location\n"
        )

        for gene in tn_coordinates:
            gene_chrom = ref_tid_roman.get(gene_coordinates.get(gene)[0])
            tncoordinates = [
                ins - chr_lengths_cumsum.get(gene_chrom)
                for ins in tn_coordinates[gene][3]
            ]

            if gene in aliases_designation:
                gene_alias = aliases_designation.get(gene)[0]
            else:
                gene_alias = gene

            f.write(
                gene_alias
                + "\t"
                + str(tn_coordinates[gene][0])
                + "\t"
                + str(tn_coordinates[gene][1] - chr_lengths_cumsum.get(gene_chrom))
                + "\t"
                + str(tn_coordinates[gene][2] - chr_lengths_cumsum.get(gene_chrom))
                + "\t"
                + str(tncoordinates)
                + "\t"
                + str(tn_coordinates[gene][4])
                + "\n"
            )


def save_per_essential_insertions(
    filename,
    tn_coordinates,
    gene_coordinates,
    chr_lengths_cumsum,
    ref_tid_roman,
    aliases_designation,
):
    """

    
    """

    with open(filename, "w") as f:

        f.write(
            "Essential gene name\tChromosome\tStart location\tEnd location\tInsertion locations\tReads per insertion location\n"
        )

        for gene in tn_coordinates:
            gene_chrom = ref_tid_roman.get(gene_coordinates.get(gene)[0])
            tncoordinates = [
                ins - chr_lengths_cumsum.get(gene_chrom)
                for ins in tn_coordinates[gene][3]
            ]

            if gene in aliases_designation:
                gene_alias = aliases_designation.get(gene)[0]
            else:
                gene_alias = gene

            f.write(
                gene_alias
                + "\t"
                + str(tn_coordinates[gene][0])
                + "\t"
                + str(tn_coordinates[gene][1] - chr_lengths_cumsum.get(gene_chrom))
                + "\t"
                + str(tn_coordinates[gene][2] - chr_lengths_cumsum.get(gene_chrom))
                + "\t"
                + str(tncoordinates)
                + "\t"
                + str(tn_coordinates[gene][4])
                + "\n"
            )

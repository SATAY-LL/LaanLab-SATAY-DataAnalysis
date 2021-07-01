import numpy as np

def get_insertions_and_reads(coordinates, tn_coordinates, readnumb_array):
    """

    """

    tn_per_gene = {}
    reads_per_gene = {}
    tn_coordinates_per_gene = {}

    for gene in coordinates:
        xx = np.where(np.logical_and(tn_coordinates[:,1] >= coordinates.get(gene)[1], tn_coordinates[:,1] <= coordinates.get(gene)[2])) #get all insertions within range of current gene
        tn_per_gene[gene] = np.size(xx)
        reads_per_gene[gene] = sum(readnumb_array[xx]) - max(readnumb_array[xx], default=0) #REMOVE LARGEST VALUE TO REDUCE NOISE

        if np.size(xx) > 0:
            tn_coordinates_per_gene[gene] = [coordinates.get(gene)[0], coordinates.get(gene)[1], coordinates.get(gene)[2], list(tn_coordinates[xx[0][0]:xx[0][-1]+1, 1]), list(readnumb_array[xx])]
        else:
            tn_coordinates_per_gene[gene] = [coordinates.get(gene)[0], coordinates.get(gene)[1], coordinates.get(gene)[2], [], []]

    return tn_per_gene, reads_per_gene, tn_coordinates_per_gene
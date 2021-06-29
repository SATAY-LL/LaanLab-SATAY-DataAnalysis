import numpy as np

def add_chromosome_length_inserts(coordinates, ref_names, chr_lengths):
    """

    """
    ll = 0
    for ii in range(1,len(ref_names)):
        ll += chr_lengths[ref_names[ii-1]]
        aa = np.where(coordinates[:,0] == ii + 1)
        coordinates[aa,1] = coordinates[aa,1] + ll

    return coordinates


def add_chromosome_length(coordinates, chr_lengths_cumsum, ref_tid_roman):
    """

    """
    for key in coordinates:
        gene_chrom = ref_tid_roman.get(coordinates.get(key)[0])
        coordinates[key][1] = coordinates.get(key)[1] + chr_lengths_cumsum.get(gene_chrom)
        coordinates[key][2] = coordinates.get(key)[2] + chr_lengths_cumsum.get(gene_chrom)
    
    return coordinates


import os
import numpy as np


def export_as_wig(wigfile, ref_names, readnumb_array, tncoordinates_array, ref_tid):
    """

    """

    readnumbwig_array = readnumb_array

    unique_index_array = np.array([], dtype=int)  # =cc
    N_uniques_perchr_list = []
    ll = 0
    for kk in ref_names:
        # get indices for current chromosome.
        index = np.where(tncoordinates_array[:, 0] == int(ref_tid[kk] + 1))

        # get all insertion locations (in tncoordinates, all rows, column 1)
        unique_index = np.unique(tncoordinates_array[index][:, 1], return_index=True)[1]

        unique_index_array = np.append(unique_index_array, (unique_index + ll), axis=0)

        # total amount unique indices found untill current chromosome
        ll += np.count_nonzero(tncoordinates_array[:, 0] == int(ref_tid[kk] + 1))
        N_uniques_perchr_list.append(ll)

    # Collect duplicates
    duplicate_list = []  # =dd
    ll = 0
    index_last_unique_previous_chromosome = 0
    for ii in N_uniques_perchr_list:
        index_last_unique = np.where(unique_index_array <= ii)[0][-1]
        for jj in range(ll, ii):
            if (
                int(jj)
                not in unique_index_array[
                    index_last_unique_previous_chromosome:index_last_unique
                ]
            ):
                duplicate_list.append(jj)
        index_last_unique_previous_chromosome = index_last_unique
        ll = ii

    # SUM READNUMB VALUES AT INDEX IN DUPLICATE_LIST AND DUPLICATE_LIST-1
    for ii in duplicate_list:
        readnumbwig_array[ii - 1] = readnumbwig_array[ii - 1] + readnumbwig_array[ii]

    tncoordinateswig_duplicatesremoved_array = np.delete(
        tncoordinates_array, duplicate_list, axis=0
    )
    readnumbwig_duplicatesremoved_array = np.delete(
        readnumbwig_array, duplicate_list, axis=0
    )

    # Write wigfile
    with open(wigfile, "w") as f:
        base = os.path.basename(wigfile)
        filename = os.path.splitext(base)[0]
        f.write("track type=wiggle_0 ,maxheightPixels=60 name=" + filename + "\n")
        for kk in ref_names:
            f.write("VariableStep chrom=chr" + kk + "\n")

            index = np.where(
                tncoordinateswig_duplicatesremoved_array[:, 0] == int(ref_tid[kk] + 1)
            )  # get indices for current chromosome.
            for ii in index[0]:
                f.write(
                    str(tncoordinateswig_duplicatesremoved_array[ii][1])
                    + " "
                    + str(readnumbwig_duplicatesremoved_array[ii])
                    + "\n"
                )

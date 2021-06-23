import numpy as np

from . import get_chromosome_names


def save_bam_as_wig(bam, bamfile: str):
    """


    Parameters
    ----------

    """

    wigfile = bamfile + ".wig"

    ref_tid = get_chromosome_names(bam)
    ref_names = list(ref_tid.keys())

    with open(wigfile, "w") as f:
        f.write("track type=wiggle_0 ,maxheightPixels=60 name=" + filename + "\n")
        for kk in ref_names:
            f.write("VariableStep chrom=chr" + kk + "\n")

            index = np.where(
                tncoordinateswig_duplicatesremoved_array[:, 0]
                == int(ref_tid_dict[kk] + 1)
            )  # get indices for current chromosome.

            for ii in index[0]:
                f.write(
                    str(tncoordinateswig_duplicatesremoved_array[ii][1])
                    + " "
                    + str(readnumbwig_duplicatesremoved_array[ii])
                    + "\n"
                )

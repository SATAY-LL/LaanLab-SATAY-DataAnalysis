import numpy as np

def correct_read_position(flags, start_position, readlength):
    """
    Correct starting position for reads with reversed orientation    

    """

    flag0coor_array = np.where(flags == 1)  # coordinates reads 5' -> 3'
    flag16coor_array = np.where(flags == -1)  # coordinates reads 3' -> 5'

    startdirect_array = start_position[flag0coor_array]
    flagdirect_array = flags[flag0coor_array]

    startindirect_array = (
        start_position[flag16coor_array] + readlength[flag16coor_array]
    )
    flagindirect_array = flags[flag16coor_array]

    start_position_corrected = np.concatenate(
        (startdirect_array, startindirect_array), axis=0
    )
    flags_corrected = np.concatenate((flagdirect_array, flagindirect_array), axis=0)

    start_sortindices = start_position_corrected.argsort(
        kind="mergesort"
    )  # use mergesort for stable sorting
    start_position_corrected = start_position_corrected[start_sortindices]
    flags_corrected = flags_corrected[start_sortindices]

    return start_position_corrected, flags_corrected

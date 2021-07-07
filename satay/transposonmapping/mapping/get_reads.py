import numpy as np
import timeit

from satay.transposonmapping.mapping.find_chromosome_reads import find_chromosome_reads
from satay.transposonmapping.mapping.correct_read_position import correct_read_position

from satay.transposonmapping.properties import (
    get_chromosome_names,
    get_sequence_length,
    get_chromosome_reads,
)


def get_reads(bam):
    """
 This function retrieves  all reads within a specified genomic region. 
     
    Usage
    ----------
    file1,file2,file3=get_reads(bam)

    Parameters
    ----------
    bam :  The output for the function pysam.AlignmentFile(bamfile, "rb")

    Returns
    -------
    readnumb_array: reads per genomic region
    tncoordinates_array: Array of three columns where the 2nd one indicated the start position where there was a transposon
    tncoordinatescopy_array: A copy from tncoordinates_array
    """

    # Get chromosome properties
    ref_tid = get_chromosome_names(bam)
    ref_names = list(ref_tid.keys())
    chr_lengths, _ = get_sequence_length(bam)
    chr_mapped_reads = get_chromosome_reads(bam)

    #%% GET ALL READS WITHIN A SPECIFIED GENOMIC REGION
    tnnumber_dict = {}
    unique_insertions = 0

    for chr_num in ref_names:
        timer_start = timeit.default_timer()

        print("Getting reads for chromosome %s ..." % chr_num)
        chromosome = bam.fetch(chr_num, 0, chr_lengths[chr_num], until_eof=True)
        N_reads = chr_mapped_reads[chr_num][2]

        [flags, start_position, readlength] = find_chromosome_reads(chromosome, N_reads)

        start_position_corrected, flags_corrected = correct_read_position(
            flags, start_position, readlength
        )

        # CREATE ARRAY OF START POSITION AND FLAGS OF ALL READS IN GENOME
        ref_tid_kk = int(ref_tid[chr_num] + 1)
        if unique_insertions == 0:
            tncoordinates_array = np.array([])

        unique_insertions_per_read = 0  # Number of unique reads per insertion
        num_transposons = 1  # Number of unique reads in current chromosome (Number of transposons in current chromosome)
        for ii in range(1, len(start_position_corrected)):
            if (
                abs(start_position_corrected[ii] - start_position_corrected[ii - 1])
                <= 2
                and flags_corrected[ii] == flags_corrected[ii - 1]
            ):  # If two subsequent reads are within two basepairs and have the same orientation, add them together.
                unique_insertions_per_read += 1
            else:
                avg_start_pos = abs(
                    round(
                        np.mean(
                            start_position_corrected[
                                ii - unique_insertions_per_read - 1 : ii
                            ]
                        )
                    )
                )
                if tncoordinates_array.size == 0:  # include first read
                    tncoordinates_array = np.array(
                        [ref_tid_kk, int(avg_start_pos), int(flags_corrected[ii - 1])]
                    )
                    readnumb_list = [unique_insertions_per_read + 1]
                else:
                    tncoordinates_array = np.vstack(
                        (
                            tncoordinates_array,
                            [
                                ref_tid_kk,
                                int(avg_start_pos),
                                int(flags_corrected[ii - 1]),
                            ],
                        )
                    )
                    readnumb_list.append(unique_insertions_per_read + 1)
                unique_insertions_per_read = 0
                num_transposons += 1
                unique_insertions += 1

            if ii == len(start_position_corrected) - 1:  # include last read
                avg_start_pos = abs(
                    round(
                        np.mean(
                            start_position_corrected[
                                ii - unique_insertions_per_read - 1 : ii
                            ]
                        )
                    )
                )
                tncoordinates_array = np.vstack(
                    (
                        tncoordinates_array,
                        [ref_tid_kk, int(avg_start_pos), int(flags_corrected[ii - 1])],
                    )
                )
                readnumb_list.append(unique_insertions_per_read + 1)

        tnnumber_dict[chromosome] = num_transposons

        timer_end = timeit.default_timer()
        print(
            "Chromosome %s completed in %.3f seconds"
            % (chromosome, (timer_end - timer_start))
        )
        print("")

    readnumb_array = np.array(readnumb_list)

    tncoordinatescopy_array = np.array(tncoordinates_array, copy=True)

    return readnumb_array, tncoordinates_array, tncoordinatescopy_array


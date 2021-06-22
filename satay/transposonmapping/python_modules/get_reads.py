import numpy as np
import timeit

from .samflag import samflags
from .get_chromosome_properties import get_chromosome_names
from .get_chromosome_properties import get_sequence_length
from .get_chromosome_properties import get_chromosome_reads


def get_reads(bam):
    """


    Parameters
    ----------
    bam : 

    Returns
    -------

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

        # REFACTOR BELOW

        # CREATE ARRAY OF START POSITION AND FLAGS OF ALL READS IN GENOME
        ref_tid_kk = int(ref_tid[chr_num] + 1)
        if unique_insertions == 0:
            tncoordinates_array = np.array([])

        unique_insertions_per_read = 0  # Number of unique reads per insertion
        num_transposons = 1  # Number of unique reads in current chromosome (Number of transposons in current chromosome)
        for ii in range(1, len(start_position_corrected)):
            if (
                abs(start_position_corrected[ii] - start_position_corrected[ii - 1]) <= 2
                and flags_corrected[ii] == flags_corrected[ii - 1]
            ):  # If two subsequent reads are within two basepairs and have the same orientation, add them together.
                unique_insertions_per_read += 1
            else:
                avg_start_pos = abs(round(np.mean(start_position_corrected[ii - unique_insertions_per_read - 1 : ii])))
                if tncoordinates_array.size == 0:  # include first read
                    tncoordinates_array = np.array(
                        [ref_tid_kk, int(avg_start_pos), int(flags_corrected[ii - 1])]
                    )
                    readnumb_list = [unique_insertions_per_read + 1]
                else:
                    tncoordinates_array = np.vstack(
                        (
                            tncoordinates_array,
                            [ref_tid_kk, int(avg_start_pos), int(flags_corrected[ii - 1])],
                        )
                    )
                    readnumb_list.append(unique_insertions_per_read + 1)
                unique_insertions_per_read = 0
                num_transposons += 1
                unique_insertions += 1

            if ii == len(start_position_corrected) - 1:  # include last read
                avg_start_pos = abs(round(np.mean(start_position_corrected[ii - unique_insertions_per_read - 1 : ii])))
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
    del readnumb_list

    tncoordinatescopy_array = np.array(tncoordinates_array, copy=True)

    return readnumb_array, tncoordinates_array, tncoordinatescopy_array


def find_chromosome_reads(chromosome, N_reads: int):
    """
    
    Syntax: [flags, start_position, readlength] = get_chromosome_reads(chromosome, N_reads)


    """
    # Initialize arrays
    start_position = np.empty(shape=(N_reads), dtype=int)
    flags = np.empty(shape=(N_reads), dtype=int)
    readlength = np.empty(shape=(N_reads), dtype=int)

    for index, reads in enumerate(chromosome):
        read = str(reads).split("\t")

        start_position[index] = int(read[3]) + 1

        # GET FLAG FOR EACH READ. IF READ ON FORWARD STRAND, ASSIGN VALUE 1, IF READ ON REVERSE STRAND ASSIGN VALUE -1, IF READ UNMAPPED OR SECONDARY ALIGNMENT ASSIGN VALUE 0
        #            flag_array[read_counter] = int(read[1])
        samprop = samflags(flag=int(read[1]), verbose=False)[1]
        if "read reverse strand" in samprop:
            flags[index] = -1
        else:
            flags[index] = 1
        if "not primary alignment" in samprop or "read unmapped" in samprop:
            flags[index] = 0

        cigarmatch_list = []
        if not reads.cigartuples == None:
            for cigar_type, cigar_length in reads.cigartuples:
                if cigar_type == 0:
                    cigarmatch_list.append(cigar_length)
                elif cigar_type == 2:
                    cigarmatch_list.append(cigar_length)
        match_length = sum(cigarmatch_list)

        readlength[index] = match_length  # int(len(read[9]))

    return flags, start_position, readlength


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

    start_position_corrected = np.concatenate((startdirect_array, startindirect_array), axis=0)
    flags_corrected = np.concatenate((flagdirect_array, flagindirect_array), axis=0)

    start_sortindices = start_position_corrected.argsort(
        kind="mergesort"
    )  # use mergesort for stable sorting
    start_position_corrected = start_position_corrected[start_sortindices]
    flags_corrected = flags_corrected[start_sortindices]

    return start_position_corrected, flags_corrected

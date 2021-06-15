import numpy as np
import timeit

from .samflag import samflags


def get_reads(bam, ref_name: list, chr_mapped_reads: dict):
    """


    Parameters
    ----------
    bam : 
    ref_name : list
    chr_mapped_reads : dict

    Returns
    -------

    """

    #%% GET ALL READS WITHIN A SPECIFIED GENOMIC REGION
    tnnumber_dict = {}
    unique_insertions = 0

    for chr_num in ref_name:
        timer_start = timeit.default_timer()

        print("Getting reads for chromosome %s ..." % chr_num)
        chromosome = bam.fetch(chr_num, 0, chr_length_dict[chr_num], until_eof=True)
        N_reads = chr_mapped_reads[chr_num][2]

        [flags, start_position, readlength] = get_chromosome_reads(chromosome, N_reads)

        start2_array, flag2_array = correct_read_position(
            flags, start_position, readlength
        )

        # REFACTOR BELOW

        # CREATE ARRAY OF START POSITION AND FLAGS OF ALL READS IN GENOME
        ref_tid_kk = int(ref_tid_dict[chr_num] + 1)
        if unique_insertions == 0:
            tncoordinates_array = np.array([])

        mm = 0  # Number of unique reads per insertion
        jj = 1  # Number of unique reads in current chromosome (Number of transposons in current chromosome)
        for ii in range(1, len(start2_array)):
            if (
                abs(start2_array[ii] - start2_array[ii - 1]) <= 2
                and flag2_array[ii] == flag2_array[ii - 1]
            ):  # If two subsequent reads are within two basepairs and have the same orientation, add them together.
                mm += 1
            else:
                avg_start_pos = abs(round(np.mean(start2_array[ii - mm - 1 : ii])))
                if tncoordinates_array.size == 0:  # include first read
                    tncoordinates_array = np.array(
                        [ref_tid_kk, int(avg_start_pos), int(flag2_array[ii - 1])]
                    )
                    readnumb_list = [mm + 1]
                else:
                    tncoordinates_array = np.vstack(
                        (
                            tncoordinates_array,
                            [ref_tid_kk, int(avg_start_pos), int(flag2_array[ii - 1])],
                        )
                    )
                    readnumb_list.append(mm + 1)
                mm = 0
                jj += 1
                unique_insertions += 1

            if ii == len(start2_array) - 1:  # include last read
                avg_start_pos = abs(round(np.mean(start2_array[ii - mm - 1 : ii])))
                tncoordinates_array = np.vstack(
                    (
                        tncoordinates_array,
                        [ref_tid_kk, int(avg_start_pos), int(flag2_array[ii - 1])],
                    )
                )
                readnumb_list.append(mm + 1)

        tnnumber_dict[chromosome] = jj

        del (
            jj,
            start_array,
            flag_array,
            readlength_array,
            flag2_array,
            start2_array,
            ref_tid_kk,
        )

        timer_end = timeit.default_timer()
        print(
            "Chromosome %s completed in %.3f seconds"
            % (chromosome, (timer_end - timer_start))
        )
        print("")

    readnumb_array = np.array(readnumb_list)
    del readnumb_list

    tncoordinatescopy_array = np.array(tncoordinates_array, copy=True)

    return readnumb_array, tncoordinatescopy_array


def get_chromosome_reads(chromosome, N_reads: int):
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

    start2_array = np.concatenate((startdirect_array, startindirect_array), axis=0)
    flag2_array = np.concatenate((flagdirect_array, flagindirect_array), axis=0)

    start2_sortindices = start2_array.argsort(
        kind="mergesort"
    )  # use mergesort for stable sorting
    start2_array = start2_array[start2_sortindices]
    flag2_array = flag2_array[start2_sortindices]

    return start2_array, flag2_array

import numpy as np

from .samflag import samflags

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

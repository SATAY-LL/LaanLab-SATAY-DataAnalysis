import warnings


def get_chromosome_names(bam):
    """

    Parameters
    ----------
    bam : dict        

    Returns
    -------
    ref_tid : dict
    
    """

    # ref_tid = {str(name): int(bam.get_tid(name)) + 1 for name in bam.get_reference_name}

    ref_tid = {} # 'I' | 0, 'II' | 1, ...
    for i in range(bam.nreferences): #if bam.nreferences does not work, use range(17) #16 chromosomes and the mitochondrial chromosome
        ref_name = bam.get_reference_name(i)
        ref_tid[ref_name] = bam.get_tid(ref_name)

    return ref_tid


def get_sequence_length(bam):
    """
    
    Parameters
    ----------
    bam : dict        

    Returns
    -------
    chr_lengths : dict
    chr_lengths_cumsum : dict

    """

    chr_lengths = {}  # 'I' | 230218, 'II' | 813184, ...
    chr_lengths_cumsum = {}  # 'I' | 0, 'II' | 230218, 'III' |  1043402, ...
    ref_summedlength = 0
    ref_tid = get_chromosome_names(bam)
    for key in ref_tid.keys():
        ref_length = bam.get_reference_length(key)
        chr_lengths[key] = bam.get_reference_length(key)
        chr_lengths_cumsum[key] = ref_summedlength
        ref_summedlength += ref_length

    return chr_lengths, chr_lengths_cumsum


def get_chromosome_reads(bam):
    """
    
    Parameters
    ----------
    bam : dict        

    Returns
    -------
    mapped_reads : dict[str, list]
        Syntax 'I' | [mapped, unmapped, total reads]

    """
    stats = bam.get_index_statistics()
    mapped_reads = {}
    for stat in stats:
        mapped_reads[stat[0]] = [stat[1], stat[2], stat[3]]
        if stat[2] != 0:
            warnings.warn("Unmapped reads found in chromosome " + stat[0])

    return mapped_reads

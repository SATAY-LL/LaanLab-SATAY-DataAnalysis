# def get_chromosome_names(bam):

#     ref_tid_dict = {} # 'I' | 0, 'II' | 1, ...
#     ref_name_list = [] # 'I', 'II', ...
#     for ref_name in bam.get_reference_name:
#         ref_tid_dict[ref_name] = bam.get_tid(ref_name)
#         ref_name_list.append(ref_name)

#     return

from typing import Dict
import warnings


def get_chromosome_names(bam) -> Dict:
    """

    Parameters
    ----------
    bam : dict        

    Returns
    -------
    ref_tid : dict
    
    """

    ref_tid = {str(name): int(bam.get_tid(name)) + 1 for name in bam.get_reference_name}

    return ref_tid


def get_sequence_length(bam):
    """
    
    Parameters
    ----------
    bam : dict        

    Returns
    -------

    """

    seq_lengths = {}  # 'I' | 230218, 'II' | 813184, ...
    cumsum_seq_lengths = {}  # 'I' | 0, 'II' | 230218, 'III' |  1043402, ...
    ref_summedlength = 0
    ref_tid = get_chromosome_names(bam)
    for key in ref_tid:
        ref_length = bam.get_reference_length(key)
        seq_lengths[key] = ref_length
        cumsum_seq_lengths[key] = ref_summedlength
        ref_summedlength += ref_length

    return seq_lengths, cumsum_seq_lengths


def get_chromosome_reads(bam) -> Dict:
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

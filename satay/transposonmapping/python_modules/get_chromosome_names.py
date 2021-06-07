# def get_chromosome_names(bam):

#     ref_tid_dict = {} # 'I' | 0, 'II' | 1, ...
#     ref_name_list = [] # 'I', 'II', ...
#     for ref_name in bam.get_reference_name:
#         ref_tid_dict[ref_name] = bam.get_tid(ref_name)
#         ref_name_list.append(ref_name)

#     return

from typing import Dict


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


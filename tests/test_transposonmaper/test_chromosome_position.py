# -*- coding: utf-8 -*-

import os
import pkg_resources

from satay.transposonmapping.python_modules.chromosome_and_gene_positions import (
    chromosome_position,
)


def test_packaged_data():
    """Test default arguments and outputs type and content """

    data_path = pkg_resources.resource_filename("satay", "data_files/")
    data_file = "Saccharomyces_cerevisiae.R64-1-1.99.gff3"
    filepath = os.path.join(data_path, data_file)

    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file=filepath)

    assert isinstance(chr_length_dict, dict), "Expected dict type"
    assert isinstance(chr_start_pos_dict, dict), "Expected dict type"
    assert isinstance(chr_end_pos_dict,dict), "Expected dict type"
    assert len(chr_start_pos_dict) == 17, "Expected dictionary of size 17"
    assert len(chr_end_pos_dict) == 17, "Expected dictionary of size 17"
    assert len(chr_length_dict) == 17, "Expected dictionary of size 17"
    assert  chr_length_dict['II']==813184, "Expected value to be present"
    assert chr_start_pos_dict['II']==230219 ,"Expected value to be present"
    assert chr_end_pos_dict[ 'II']==1043402, "Expected value to be present"
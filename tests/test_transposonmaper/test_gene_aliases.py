# -*- coding: utf-8 -*-

import os
import pkg_resources

from satay.transposonmapping.python_modules.gene_names import (
    gene_aliases,
)


def test_packaged_data():
    """Test default arguments and outputs type and content """

    data_path = pkg_resources.resource_filename("satay", "data_files/")
    data_file = 'Yeast_Protein_Names.protein_names.txt'
    filepath = os.path.join(data_path, data_file)

    aliases_designation_dict, aliases_sgd_dict, aliases_swissprot_dict = gene_aliases(gene_information_file=filepath)

    assert isinstance(aliases_designation_dict, dict), "Expected dict type"
    assert isinstance(aliases_sgd_dict, dict), "Expected dict type"
    assert isinstance(aliases_swissprot_dict,dict), "Expected dict type"
    assert len(aliases_designation_dict) == 6725, "Expected dictionary of size 6725"
    assert len(aliases_sgd_dict) == 6725, "Expected dictionary of size 6725"
    assert len(aliases_swissprot_dict) == 6725, "Expected dictionary of size 6725"
    assert 'YMR056C' in aliases_designation_dict, "Expected YBR247C key to be present"
    assert aliases_designation_dict['YMR056C']==['AAC1'], "Expected value to be present"
# -*- coding: utf-8 -*-


from satay.transposonmapping.python_modules.chromosome_and_gene_positions import (
    chromosomename_roman_to_arabic,
)


def test_output_data():
    """Test default output """

    
    dict_1,dict_2 = chromosomename_roman_to_arabic()

    assert isinstance(dict_1, dict), "Expected dict type"
    assert isinstance(dict_2, dict), "Expected dict type"
    assert len(dict_1) ==17, "Expected dictionary of size 17"
    assert len(dict_2) ==17, "Expected dictionary of size 17"
    assert 1 in dict_1, "Expected 1 key to be present"
    assert 'I' in dict_1[1], "Expected I as translation to 1"
    assert 'I' in dict_2, "Expected I key to be present"
    assert  dict_2['I']==1, "Expected 1 as translation to I"

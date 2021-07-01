import os
import pkg_resources
import pytest
import numpy as np

import pysam
from satay.transposonmapping.properties import get_sequence_length


@pytest.fixture
def bam():
    """
    Load bamfile for testing
    """
    data_path = pkg_resources.resource_filename("satay", "data_files/files4test/")
    filename = "SRR062634.filt_trimmed.sorted.bam"
    bamfile = os.path.join(data_path, filename)
    bam = pysam.AlignmentFile(bamfile, "rb")  # open bam formatted file for reading
    return bam


def test_output_is_dict(bam):
    [sum_dict, cumsum_dict] = get_sequence_length(bam)
    assert isinstance(sum_dict, dict), "Expected dict type"
    assert isinstance(cumsum_dict, dict), "Expected dict type"


def test_output_size(bam):
    [sum_dict, cumsum_dict] = get_sequence_length(bam)
    assert len(sum_dict) == 17, "Expected dictionary of size 17"
    assert len(cumsum_dict) == 17, "Expected dictionary of size 17"


def test_output_values_cumsum(bam):
    [sum_dict, cumsum_dict] = get_sequence_length(bam)
    array = np.cumsum(list(sum_dict.values()))
    exp_values = np.hstack((np.zeros(1, dtype=int), array[:-1])).tolist()
    act_values = list(cumsum_dict.values())
    assert act_values == exp_values, "Expected cumulative sum"


def test_output_dict_keys(bam):
    [sum_dict, cumsum_dict] = get_sequence_length(bam)
    exp_key = "ref|NC_001133|"

    act_key = list(sum_dict.keys())[0]
    assert exp_key == act_key, "Expected different reference key"

    act_key = list(cumsum_dict.keys())[0]
    assert exp_key == act_key, "Expected different reference key"

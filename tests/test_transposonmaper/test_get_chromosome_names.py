import os
import pkg_resources
import pytest
import numpy as np

import pysam
from satay.transposonmapping.properties.get_chromosome_properties import (
    get_chromosome_names,
)


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


def test_output_type_is_dict(bam):
    ref_tid = get_chromosome_names(bam)
    assert isinstance(ref_tid, dict), "Expected dict type"


def test_output_dict_size(bam):
    ref_tid = get_chromosome_names(bam)
    assert len(ref_tid) == 17, "Expected dictionary of size 17"


def test_output_dict_values(bam):
    ref_tid = get_chromosome_names(bam)
    exp_values = np.arange(17).tolist()
    act_values = list(ref_tid.values())
    assert exp_values == act_values, "Expected values to range from 0 to 16"


def test_output_dict_keys(bam):
    ref_tid = get_chromosome_names(bam)
    exp_key = "ref|NC_001133|"
    act_key = list(ref_tid.keys())[0]
    assert exp_key == act_key, "Expected different reference key"

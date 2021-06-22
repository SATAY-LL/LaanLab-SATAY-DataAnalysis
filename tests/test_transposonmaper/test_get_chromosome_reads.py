import os
import pkg_resources
import pytest
import numpy as np

import pysam
from satay.transposonmapping.python_modules import get_chromosome_reads

@pytest.fixture
def bam():
    """
    Load bamfile for testing
    """
    data_path = pkg_resources.resource_filename("satay", "data_files/files4test/")
    filename = 'SRR062634.filt_trimmed.sorted.bam'       
    bamfile = os.path.join(data_path, filename)
    bam = pysam.AlignmentFile(bamfile, 'rb') #open bam formatted file for reading
    return bam 

def test_output_is_dict(bam):
    reads = get_chromosome_reads(bam)
    assert isinstance(reads, dict), "Expected dict type"

def test_output_size(bam):
    reads = get_chromosome_reads(bam)
    assert len(reads) == 17, "Expected dictionary of size 17"

def test_output_dict_keys(bam):
    reads = get_chromosome_reads(bam)
    exp_key = "ref|NC_001133|"
    
    act_key = list(reads.keys())[0]
    assert exp_key == act_key, "Expected different reference key"

def test_output_dict_values(bam):
    reads = get_chromosome_reads(bam)
    exp_value = [9, 0, 9]
    act_value = list(reads.values())[0]
    assert act_value == exp_value, "Expected different values for first key"
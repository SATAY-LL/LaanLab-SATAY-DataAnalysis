import os
import pkg_resources
import pysam
import numpy 

from satay.transposonmapping.mapping import get_reads

default_path = pkg_resources.resource_filename("satay", "data_files/files4test")
bamfile = os.path.join(default_path,
                         "SRR062634.filt_trimmed.sorted.bam")
        
bam = pysam.AlignmentFile(bamfile, "rb")

readnumb_array, tncoordinates_array, tncoordinatescopy_array=get_reads(bam)

def test_ouput_format():
    assert type(readnumb_array)==numpy.ndarray, "The outputs of get_reads.py should be  arrays"
    assert type(tncoordinates_array)==numpy.ndarray, "The outputs of get_reads.py should be  arrays"
    assert type(tncoordinatescopy_array)==numpy.ndarray, "The outputs of get_reads.py should be  arrays"

    assert tncoordinates_array.shape[1]==3, "This array should have three columns"
    assert tncoordinatescopy_array.shape[1]==3, "This array should have three columns"
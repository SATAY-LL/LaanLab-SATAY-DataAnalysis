import os
import pkg_resources
import glob

from satay.transposonmapping.transposonmapping_satay import (
    transposonmapper 
)
def test_transposonmapper_output_data():
    """Test default number of outputs"""

    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
   
    filepath = os.path.join(data_path, filename)

    transposonmapper(bamfile=filepath)
    data=os.listdir(data_path)
    
    assert len(data)>=2+6, "Expected 6 extra data files"
    
def test_transposonmapper_type_output():
    
    """ Test default type of outputs """
    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
   
    filepath = os.path.join(data_path, filename)

    transposonmapper(bamfile=filepath)
    files_txt=glob.glob(data_path + '/**/*'+ '.txt', recursive=True)
    file_wig=glob.glob(data_path + '/**/*'+ '.wig', recursive=True)
    file_bed=glob.glob(data_path + '/**/*'+ '.bed', recursive=True)
    
    assert len(files_txt)>=4 , "We need 4 txt files as output per datafile"
    assert len(file_wig)>=1 , "We need one wig file as output per datafile"
    assert len(file_bed)>=1, "We need one bed file as output per datafile"
    
    
   
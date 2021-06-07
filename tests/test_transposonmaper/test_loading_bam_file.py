
import os
from satay.transposonmapping.python_modules.loading_bam_file import loading_bam_file

def test_loading_bam_file():
    
    path = os.path.join('satay/data_files/data_merged_wt/')
    filename = 'WT_merged-techrep-a_techrep-b_trimmed.sorted.bam'
    bamfile = os.path.join(path,filename)

    output=loading_bam_file(file_path=None)
     
    assert output==bamfile, "Should be %s" % bamfile

    assert loading_bam_file(file_path)==file_path , "Should be %s" % file_path

    
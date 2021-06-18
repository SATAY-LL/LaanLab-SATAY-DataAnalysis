import os
import pkg_resources
import glob
import pandas as pd

from satay.transposonmapping.transposonmapping_satay import (
    transposonmapper 
)


    
def test_output_format_pergene():
    """Test default number of outputs"""

    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
       
    filepath = os.path.join(data_path, filename)
    
    transposonmapper(bamfile=filepath)
    
    data_list=os.listdir(data_path)
    pergene_file=glob.glob(data_path + '/**/*'+ '_pergene.txt', recursive=True)
    
    data=pd.read_csv(pergene_file[0],delimiter='\t')
    
    
    assert data.columns.tolist() == ['Gene name', 
                                     'Number of transposons per gene',
                                     'Number of reads per gene'] , "wrong columns names"
    assert len(data)==6600, "Incorrect length of the file"
    assert data.iloc[:,0].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,1].dtypes=='int64', "Incorrect data type"
    assert data.iloc[:,2].dtypes=='int64', "Incorrect data type"
    
    
def test_output_format_peressential_insertions():
    """Test default number of outputs"""
    
    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
    filepath = os.path.join(data_path, filename)
    transposonmapper(bamfile=filepath)
    
    data_list=os.listdir(data_path)
    
    pergene_insert_file=glob.glob(data_path + '/**/*'+ '_pergene_insertions.txt', recursive=True)
    
    data=pd.read_csv(pergene_insert_file[0],delimiter='\t')
    
    
    assert data.columns.tolist() == ['Gene name',
                                     'Chromosome',
                                     'Start location',
                                     'End location',
                                     'Insertion locations',
                                     'Reads per insertion location'] , "wrong columns names"
    assert len(data)==6600, "Incorrect length of the file"
    assert data.iloc[:,0].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,1].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,2].dtypes=='int64', "Incorrect data type"
    assert data.iloc[:,3].dtypes=='int64', "Incorrect data type"
    assert data.iloc[:,4].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,5].dtypes=='object' , "Incorrect data type"
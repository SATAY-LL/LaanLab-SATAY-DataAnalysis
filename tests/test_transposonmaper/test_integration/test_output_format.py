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
    

    file=glob.glob(data_path + '/**/*'+ '_pergene.txt', recursive=True)
    
    data=pd.read_csv(file[0],delimiter='\t')
    
    
    assert data.columns.tolist() == ['Gene name', 
                                     'Number of transposons per gene',
                                     'Number of reads per gene'] , "wrong columns names"
    assert len(data)==6600, "Incorrect length of the file"
    assert data.iloc[:,0].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,1].dtypes=='int64', "Incorrect data type"
    assert data.iloc[:,2].dtypes=='int64', "Incorrect data type"
    
    
def test_output_format_pergene_insertions():
    """Test default number of outputs"""
    
    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
    filepath = os.path.join(data_path, filename)
    transposonmapper(bamfile=filepath)
    
        
    file=glob.glob(data_path + '/**/*'+ '_pergene_insertions.txt', recursive=True)
    
    data=pd.read_csv(file[0],delimiter='\t')
    
    
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
    
    
def test_output_format_peressentials():
    """Test default number of outputs"""

    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
       
    filepath = os.path.join(data_path, filename)
    
    transposonmapper(bamfile=filepath)
    
   
    file=glob.glob(data_path + '/**/*'+ '_peressential.txt', recursive=True)
    
    data=pd.read_csv(file[0],delimiter='\t')
    
    
    assert data.columns.tolist() == ['Gene name', 
                                     'Number of transposons per gene',
                                     'Number of reads per gene'] , "wrong columns names"
    #assert len(data)==1186, "Incorrect length of the file"
    assert data.iloc[:,0].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,1].dtypes=='int64', "Incorrect data type"
    assert data.iloc[:,2].dtypes=='int64', "Incorrect data type"
    
def test_output_format_peressentials_insertions():
    """Test default number of outputs"""
    
    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
    filepath = os.path.join(data_path, filename)
    transposonmapper(bamfile=filepath)
    
        
    file=glob.glob(data_path + '/**/*'+ '_peressential_insertions.txt', recursive=True)
    
    data=pd.read_csv(file[0],delimiter='\t')
    
    
    assert data.columns.tolist() == ['Essential gene name',
                                     'Chromosome',
                                     'Start location',
                                     'End location',
                                     'Insertion locations',
                                     'Reads per insertion location'] , "wrong columns names"
    #assert len(data)==1186, "Incorrect length of the file"
    assert data.iloc[:,0].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,1].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,2].dtypes=='int64', "Incorrect data type"
    assert data.iloc[:,3].dtypes=='int64', "Incorrect data type"
    assert data.iloc[:,4].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,5].dtypes=='object' , "Incorrect data type"
    
def test_output_format_wig():
    
    """Test default number of outputs"""
    
    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
    filepath = os.path.join(data_path, filename)
    transposonmapper(bamfile=filepath)
    
        
    file=glob.glob(data_path + '/**/*'+ '.wig', recursive=True)
    
    data=pd.read_csv(file[0],delimiter= ' ')
    
    assert data.columns.tolist() == ['track',
                                     'type=wiggle_0',
                                     ',maxheightPixels=60',
                                     'name=SRR062634.filt_trimmed.sorted.bam'] , "wrong columns names"
    assert data.iloc[:,0].dtypes=='object' , "Incorrect data type"
    assert data.iloc[:,1].dtypes=='object' , "Incorrect data type"
    assert data.iloc[0,1] =='chrom=chrref|NC_001133|' , "Incorrect first headline"
    
    
def test_output_format_bed():
    
    """Test default number of outputs"""
    
    data_path = pkg_resources.resource_filename("satay", "data_files/files4test")
    filename = 'SRR062634.filt_trimmed.sorted.bam'
    filepath = os.path.join(data_path, filename)
    transposonmapper(bamfile=filepath)
    
        
    file=glob.glob(data_path + '/**/*'+ '.bed', recursive=True)
    
    data=pd.read_csv(file[0],delimiter= ' ')
    
    assert data.columns.tolist() == ['track', 
                                     'name=SRR062634.filt_trimmed.sorted.bam', 
                                     'useScore=1'], "wrong columns names"
    assert data.iloc[:,0].dtypes=='int64' , "Incorrect data type"
    assert data.iloc[:,1].dtypes=='object' , "Incorrect data type"
    assert data.index[0][0] =='chrref|NC_001133|' or data.index[0][0] =='chrI' , "Incorrect first headline"
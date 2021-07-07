
from satay.transposonmapping.importing import (
    load_default_files,read_genes
)

def test_output_format():
    a,b,c=load_default_files(gff_file=None,essentials_file=None,gene_names_file=None)
    a_0,b_0,c_0=read_genes(gff_file=a,essentials_file=b,gene_names_file=c)
    
    assert type(a_0)==dict, "the gene coordinates have to be a dict"
    assert type(b_0)==dict, "the gene coordinates have to be a dict"
    assert type(c_0)==dict, "the gene coordinates have to be a dict"

def test_output_length():
    a,b,c=load_default_files(gff_file=None,essentials_file=None,gene_names_file=None)
    a_0,b_0,c_0=read_genes(gff_file=a,essentials_file=b,gene_names_file=c)
    
    assert len(a_0)>=6600, "the total number of genes should not be less than 6600"
    assert len(b_0)<6600, "the total number of essential genes should not be more than the number of genes"
    assert len(c_0)>=6600, "the total number of genes should not be less than 6600"
    
    

def test_output_content_gff():
    a,b,c=load_default_files(gff_file=None,essentials_file=None,gene_names_file=None)
    a_0,b_0,c_0=read_genes(gff_file=a,essentials_file=b,gene_names_file=c)
    
    #read the first value of the dict
    first_value=next(iter(a_0.values()))
    
    # read the first key
    first_key=next(iter(a_0))
    
    assert first_value==['I', 335, 649, '+'], "The first value of the gene coordinates is wrong"
    
    assert first_key== 'YAL069W', "The first gene in the array should be YAL069W"
    
    
def test_output_content_essentials():
    a,b,c=load_default_files(gff_file=None,essentials_file=None,gene_names_file=None)
    a_0,b_0,c_0=read_genes(gff_file=a,essentials_file=b,gene_names_file=c)
    
    #read the first value of the dict
    first_value=next(iter(b_0.values()))
    
    # read the first key
    first_key=next(iter(b_0))
    
    assert first_value==['I', 147594, 151166, '-'], "The first value of the gene coordinates is wrong"
    
    assert first_key== 'YAL001C', "The first gene in the array should be YAL001C"
    
def test_output_content_names():
    a,b,c=load_default_files(gff_file=None,essentials_file=None,gene_names_file=None)
    a_0,b_0,c_0=read_genes(gff_file=a,essentials_file=b,gene_names_file=c)
    
    #read the first value of the dict
    first_value=next(iter(c_0.values()))
    
    # read the first key
    first_key=next(iter(c_0))
    
    assert first_value==['AAC1'], "The first value of the gene names is wrong"
    
    assert first_key== 'YMR056C', "The first gene in the array should be YMR056C"
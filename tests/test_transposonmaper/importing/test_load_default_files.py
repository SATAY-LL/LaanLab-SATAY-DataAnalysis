import os
import pkg_resources
import glob
import pandas as pd

from satay.transposonmapping.importing import (
    load_default_files,load_sgd_tab
)

def test_load_default_files():
   inputs=["a.gff3","b.txt","c.txt"]
   a, b,c= load_default_files(inputs[0],inputs[1],inputs[2])
   assert a==inputs[0] , "expect output to be the same as input is not None"
   assert b==inputs[1] , "expect output to be the same as input is not None"
   assert c==inputs[2] , "expect output to be the same as input is not None"

def test_load_default_files_None():
    a,b,c= load_default_files(gff_file=None,essentials_file=None,gene_names_file=None)
    default_path = pkg_resources.resource_filename("satay", "data_files/")

    gff_file = os.path.join(
            default_path, "Saccharomyces_cerevisiae.R64-1-1.99.gff3"
        )
    essentials_file = os.path.join(
            default_path, "Cerevisiae_AllEssentialGenes_List.txt"

        )
    gene_names_file = os.path.join(default_path, "Yeast_Protein_Names.txt")

    assert a==gff_file , "it should be this file if the input is None"
    assert b==essentials_file , "it should be this file if the input is None"
    assert c==gene_names_file , "it should be this file if the input is None"


def test_load_sgd_tab():
    inputs=["a.tab"]
    a=load_sgd_tab(inputs[0])
    assert a==inputs[0] , "expect output to be the same as input is not None"

def test_load_sgd_tab_None():
    a=load_sgd_tab(sgd_features_file=None)
    default_path = pkg_resources.resource_filename("satay", "data_files/")

    sgd_features_file = os.path.join(
            default_path, "SGD_features.tab"
        )
    assert a==sgd_features_file, "it should be this file if the input is None"
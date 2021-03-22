# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:39:53 2021

@author: gregoryvanbeek

This script removes transposon insertions in .bed and .wig files that were mapped outside the chromosomes, creates consistent naming for chromosomes and change the header of files with custom headers.
"""

import os, sys

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import chromosome_position
from chromosome_names_in_files import chromosome_name_bedfile, chromosome_name_wigfile



#%% INPUT -> FILEPATH TAKES PATHS TO .BED OR .WIG FILES
filepath = r""
custom_header = ""
split_chromosomes = False



#%%
def cleanfiles(filepath=None, custom_header=None, split_chromosomes=False):
    '''
    This code reads a .bed or .wig file and remove any insertions that were mapped outside a chromosome.
    Mapping of a read outside a chromosome can happen during the alignment and transposon mapping steps and means that the position of an insertions site of a read is larger than the length of the chromosome it is mapped to.
    This function creates a new file with the same name as the inputfile with the extension _clean.bed or _clean.wig.
    This is saved at the same location as the input file.
    In this _clean file the redundant insertions that were mapped outside the chromosome are removed.
    The lengths of the chromosomes are determined the python function 'chromosome_position' which is part of the python module 'chromosome_and_gene_positions.py'.
    This module gets the lengths of the chromosomes from a .gff file downloaded from SGD (https://www.yeastgenome.org/).
    Besides removing the reads outside the chromosomes, it also changes the names of the chromosomes to roman numerals and a custom header can be inputted (optional).
    Finally, the bed and wig files can be split up in separate files for each chromosome. These are placed in _chromosomesplit folder located at the location of the bed or wig file.
    '''

    if filepath == None:
        sys.exit(0)
    else:
        assert os.path.isfile(filepath), 'File not found: %s' % filepath

    chr_length_dict = chromosome_position()[0]

    filepath_splitext = os.path.splitext(filepath)
    exten = filepath_splitext[1]



    num_roman = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']



    if exten == ".bed":
        print("Bed file loaded %s" % filepath)

        chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict = chromosome_name_bedfile(filepath)

        with open(filepath, "r") as f:
            lines = f.readlines()


        with open(filepath_splitext[0]+"_clean.bed", "w") as w:
            #write header
            if custom_header == None or custom_header == "":
                w.write(lines[0])
            else:
                w.write("track name=" + str(custom_header) + " useScore=1\n")

            for chrom in num_roman:
                print("evaluating chromosome %s" % chrom)

                for line in lines[chrom_start_line_dict.get(chrom): chrom_end_line_dict.get(chrom)+1]:
                    line_list = " ".join(line.strip("\n").split()).split(" ")
                    if int(line_list[2]) > chr_length_dict.get(chrom) or int(line_list[1]) < 0:
                        print("Line removed: %s" % line)
                    else:
                        for romanname, chromname in chrom_names_dict.items():
                            if chromname == line_list[0].replace("chr",""):
                                chrom_nameroman = romanname
                        w.write("chr" + str(chrom_nameroman) + " " + str(line_list[1]) + " " + str(line_list[2]) + " " + str(line_list[3]) + " " + str(line_list[4]) + "\n")

            
            for line in lines[chrom_end_line_dict.get("XVI")+1:]:
                line_list = " ".join(line.strip("\n").split()).split(" ")
                w.write("chrM" + " " + str(line_list[1]) + " " + str(line_list[2]) + " " + str(line_list[3]) + " " + str(line_list[4]) + "\n")


        if split_chromosomes == True:
            path = os.path.dirname(filepath)
            name = os.path.splitext(os.path.basename(filepath_splitext[0]+"_clean.bed"))[0]

            directoryname = os.path.join(path, name + '_chromosomesplit')

            if not os.path.exists(directoryname):
                os.mkdir(directoryname)

            chromosome_names = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']

            chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict = chromosome_name_bedfile(os.path.join(path, name+".bed"))

            with open(os.path.join(path,name+".bed"), 'r') as f:
                lines = f.readlines()
            header = lines[0]


            for chrom in chromosome_names:
                outputfile = os.path.join(directoryname, name + '_' + str(chrom) + '.bed')
                with open(outputfile, 'w+') as f:
                    f.write(header)
                    for l in range(chrom_start_line_dict.get(chrom), chrom_end_line_dict.get(chrom)+1):
                        f.write(lines[l])

                outputfile = os.path.join(directoryname, name + '_M.bed')
                with open(outputfile, 'w+') as f:
                    f.write(header)
                    for l in range(chrom_end_line_dict.get(chromosome_names[-1])+1, len(lines)):
                        f.write(lines[l])




    elif exten == ".wig":
        print("Wig file loaded %s" % filepath)

        chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(filepath)

        with open(filepath, 'r') as f:
            lines = f.readlines()

        with open(filepath_splitext[0]+"_clean.wig", "w") as w:
            #write header
            if custom_header == None:
                w.write(lines[0].replace(',',''))
            else:
                w.write("track type=wiggle_0 maxheightPixels=60 name=" + str(custom_header) + "\n")

            for chrom in num_roman:
                print("evaluating chromosome %s" % chrom)

                #replace chromosome names from reference genome with roman numerals
                chrom_headerline = lines[chrom_start_line_dict.get(chrom) - 1]
                chrom_nameline = chrom_headerline.split("=")[1].strip("\n").replace("chr","")
                for romanname, chromname in chrom_names_dict.items():
                    if chromname.replace("chr","") == chrom_nameline:
                        chrom_nameroman = romanname
                w.write("variablestep chrom=chr" + str(chrom_nameroman) + "\n") #write header for each chromosome
                for line in lines[chrom_start_line_dict.get(chrom): chrom_end_line_dict.get(chrom)]: #no '+1' in for loop, this is only for bed file
                    line_list = " ".join(line.strip("\n").split()).split(" ")
                    if int(line_list[0]) > chr_length_dict.get(chrom) or int(line_list[0]) < 0:
                        print("Line removed: %s" % line)
                    else:
                        w.write(line)


            w.write("variablestep chrom=chrM\n")
            for line in lines[chrom_end_line_dict.get("XVI")+1:]:
                w.write(line)



        if split_chromosomes == True:
            path = os.path.dirname(filepath)
            name = os.path.splitext(os.path.basename(filepath_splitext[0]+"_clean.wig"))[0]

            directoryname = os.path.join(path, name + '_chromosomesplit')

            if not os.path.exists(directoryname):
                os.mkdir(directoryname)

            chromosome_names = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','M']

            chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(os.path.join(path, name+".wig"))

            with open(os.path.join(path,name+".wig"), 'r') as f:
                lines = f.readlines()
            header = lines[0]


            for chrom in chromosome_names:
                outputfile = os.path.join(directoryname, name + '_' + str(chrom) + '.wig')
                with open(outputfile, 'w+') as f:
                    f.write(header)
                    for l in range(chrom_start_line_dict.get(chrom)-1, chrom_end_line_dict.get(chrom)):
                        f.write(lines[l])

    else:
        print("Extension not recognized")



#%%
if __name__ == '__main__':
    cleanfiles(filepath=filepath, custom_header=custom_header, split_chromosomes=split_chromosomes)


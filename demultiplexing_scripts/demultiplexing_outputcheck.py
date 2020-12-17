# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 13:46:31 2020

@author: gregoryvanbeek
"""


import os, sys

def check_fastq_pairedinterleaved(inputfile='', barcodes=[]):

    if not os.path.isfile(inputfile):
        print('WARNING: input file not found.')
        sys.exit()
    if barcodes == []:
        print('WARNING: barcodes list is empty')
        sys.exit()


    print("Assessing barcodes: ", barcodes)


    total_line_counter = 0
    header_counter = 0
    pair_counter = 0
    barcode_in_line = False
    with open(inputfile,'r') as f:
        for line in f:
            total_line_counter += 1
            if line.startswith('@'):
                header_counter += 1

                if pair_counter == 0:
                    header1 = line.split(' ')[0]
                    pair_counter += 1
                elif pair_counter == 1:
                    header2 = line.split(' ')[0]
                    if not header1 == header2:
                        print('Header %s is not the same as the previous header %s' % (header2, header1))
                    pair_counter = 0


            if not barcodes == []:
                if line.startswith('@'):
                    nextline = next(f) #NOTE: this skips the sequence line in the next loop. 
                    total_line_counter += 1

                    for barcode in barcodes:
                        if barcode.upper() in nextline:
                            barcode_in_line = True
                            break
                    if not barcode_in_line == True:
                        print('No barcode found in read %s' % header1)
                    else:
                        barcode_in_line = False


    print(total_line_counter)
    print(header_counter)
    lines_divisible_by_four = (total_line_counter / 4).is_integer()
    if lines_divisible_by_four == False:
        print("WARNING: Number of lines in fastq file is not a multiple of 4. Format of fastq file may be corrupted.")
    
    lines_divisible_by_headercounter = total_line_counter / header_counter
    if not int(lines_divisible_by_headercounter) == 4:
        print("WARNING: Number of header in fastq file does not match the number of lines in the file. Format of fastq file may be corrupted.")


#Check if length sequence lines is the same as the length of the phred score lines.
#Check if every fourth line startswith @



if __name__ == '__main__':
#    check_fastq_pairedinterleaved(inputfile=r"/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/wt1_enzo_dataset_demultiplexed_interleaved/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_pairs.fq",
#                              barcodes=["GCCACATA", "GCGAGTAA"])
#    check_fastq_pairedinterleaved(inputfile=r"/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/wt1_enzo_dataset_demultiplexed_interleaved/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_pairs.fq",
#                              barcodes=["GAGCTGAA", "GATAGACA"])
    check_fastq_pairedinterleaved(inputfile=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\wt1_dataset_enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_pairs_400lines.fq",
                                  barcodes=["GAGCTGAA", "GATAGACA"])

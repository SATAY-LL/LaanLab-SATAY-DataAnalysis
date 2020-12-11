# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 13:46:31 2020

@author: gregoryvanbeek
"""


import os, sys

def check_headers_pairedfastq(inputfile='', barcodes=[]):

    if not os.path.isfile(inputfile):
        print('WARNING: input file not found.')
        sys.exit()
    if barcodes == []:
        print('WARNING: barcodes list is empty')
        sys.exit()


    print("Assessing barcodes: ", barcodes)


    header_counter = 0
    barcode_in_line = False
    with open(inputfile,'r') as f:
        for line in f:
            if line.startswith('@'):
                if header_counter == 0:
                    header1 = line.split(' ')[0]
                    header_counter += 1
                elif header_counter == 1:
                    header2 = line.split(' ')[0]
                    if not header1 == header2:
                        print('Header %s is not the same as the previous header %s' % (header2, header1))
                    header_counter = 0

            if not line.startswith('@') and not line.startswith('+') and not 'F' in line:
                for barcode in barcodes:
                    if barcode in line:
                        barcode_in_line = True
                        break
                if not barcode_in_line == True:
                    print('No barcode found in read %s' % header1)
                else:
                    barcode_in_line = False





if __name__ == '__main__':
    check_headers_pairedfastq(inputfile=r"C:\Users\gregoryvanbeek\Documents\test_pairedinterleaved.fq",
                              barcodes=["GCCACATA", "GCGAGTAA"])

import os

def save_as_bed(bamfile, tncoordinates_array, ref_tid, readnumb_array):
    """

    """


    bedfile = bamfile + '.bed'    
    filename = os.path.basename(bamfile)

    print('Writing bed file at: ', bedfile)
    print('')


    with open(bedfile, 'w') as f:
        
        f.write('track name=' + filename + ' useScore=1\n')
        
        coordinates_counter = 0
        for tn in tncoordinates_array:
            refname = [key for key, val in ref_tid.items() if val == tn[0] - 1][0]
            if refname == 'Mito':
                refname = 'M'
            f.write('chr' + refname + ' ' + str(tn[1]) + ' ' + str(tn[1] + 1) + ' . ' + str(100+readnumb_array[coordinates_counter]*20) + '\n')
            coordinates_counter += 1
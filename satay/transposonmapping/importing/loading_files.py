import os, sys,glob

def access_files(file_path,extension):

    # EXTENSIONs : .sorted.bam, .gff3, .protein_names.txt, .essential_genes.txt
    if file_path is None:
        path='satay/data_files/'
        files=glob.glob(path + '/**/*'+ extension, recursive=True)
        output_file=files[0]
       
                
        
    else:
        output_file = os.path.join(file_path,extension)
        
    assert os.path.isfile(output_file), 'File not found at: %s' % output_file #check if given bam file exists

    return output_file

import os, sys
def loading_bam_file(file_path):
    
    if file_path is None:
        path = os.path.join('satay/data_files/data_merged_wt/')
        # filename = 'E-MTAB-4885.WT2.bam'
        filename = 'WT_merged-techrep-a_techrep-b_trimmed.sorted.bam'
        bamfile = os.path.join(path,filename)
    else:
        filename = os.path.basename(file_path)
        path = bamfile.replace(filename,'')

## - [] test whether the file has .sorted.bam extension
## - [] test whether the file exist in the provided location
    assert os.path.isfile(bamfile), 'Bam file not found at: %s' % bamfile #check if given bam file exists
   
    return bamfile
import pandas as pd
import numpy as np
from input_reads_output_fitness import fitness

genes_start_stop = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
fitness_insertions = fitness(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\Processed_data_WT1.wig')


genes_start_stop['fitness'] = ''
arr = np.zeros((len(genes_start_stop)), dtype=object)
genes_start_stop['fitness'] = arr

genes_start_stop.at[1,'fitness'].append(10)

# for i in range(0,len(fitness_insertions)):
#     for j in range (0,len(genes_start_stop)):
#         if fitness_insertions.at[i,'chromosome']== genes_start_stop.at[j,'chromosome'] and fitness_insertions.at[i,'tn start position'] <= genes_start_stop.at[j,'stop bp'] and fitness_insertions.at[i,'tn start position'] >= genes_start_stop.at[j,'start bp'] :

#             a = fitness_insertions.at[i,'fitness']
#             np.append(genes_start_stop.at[j,'fitness'],a)

#             break
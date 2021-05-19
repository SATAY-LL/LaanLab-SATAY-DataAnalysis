import pandas as pd
from input_reads_output_fitness import fitness

def fitness_eachgene(filepath_and_name):
    genes_start_stop = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
    fitness_insertions = fitness(r'filepath_and_name')
    
    genes_start_stop['fitness'] = ''
    genes_start_stop['fitness']= genes_start_stop['fitness'].astype('object')
    
    for i in range (0, len(genes_start_stop)):
        genes_start_stop.at[i,'fitness'] = []
    
    genes_start_stop = genes_start_stop.sort_values('chromosome')
    genes_start_stop = genes_start_stop.reset_index()
    for i in range(0,len(fitness_insertions)):
        print(i)
        if fitness_insertions.at[i,'chromosome'] == 'chrI':
            for j in range (0, 125):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrII':
            for j in range (125, 606):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrIII':
            for j in range (607, 809):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrIV':
            for j in range (810, 1697):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrIX':
            for j in range (1698, 1958):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrV':
            for j in range (1959, 2316):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrVI':
            for j in range (2317, 2471):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrVII':
            for j in range (2472, 3111):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrVIII':
            for j in range (3112, 3454):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrX':
            for j in range (3455, 3889):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrXI':
            for j in range (3890, 4259):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrXII':
            for j in range (4260, 4899):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrXIII':
            for j in range (4900, 5450):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrXIV':
            for j in range (5451, 5909):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrXV':
            for j in range (5910, 6548):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrXVI':
            for j in range (6549, 7093):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
        if fitness_insertions.at[i,'chromosome'] =='chrmt':
            for j in range (7094, 7148):
                if fitness_insertions.at[i, 'tn start position'] <= genes_start_stop.at[j, 'stop bp'] and fitness_insertions.at[i, 'tn start position'] >= genes_start_stop.at[j, 'start bp']:
                    a = fitness_insertions.at[i, 'fitness']
                    genes_start_stop.at[j, 'fitness'].append(a)
                    break
                
    genes_start_stop = genes_start_stop.sort_values('index')
    genes_start_stop.reset_index()
    return(genes_start_stop)
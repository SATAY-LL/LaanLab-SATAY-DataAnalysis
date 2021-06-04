In this folder are the scripts, modules and datasets needed to compute fitness values from the reads and insertions from SATAY data , in this case for WT and dnrp1 strain.
Each script generates an excel file that is needed for the next one , the order is as follows:
The dataset for the 1st script is the outcome of the python script ../genomicfeatures_dataframe.py per strain. 
1. script_basic_data_analysis_satay_libraries.py this script requires an available dataset inside the datasets folder
2. script_fitness-from-reads-per-transposon-satay-data.py generates the fitness values from the intergenic model that requires an excel file from the previous script
3. script_checking_interactors_from_intergenic_fitness_model.py which checks the estimated interactors of nrp1 with existing constanzo interactors. 

**It is still no tested code**

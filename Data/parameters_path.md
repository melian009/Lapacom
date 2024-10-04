# List of file path which contains the parameter used in the simulation on life cycles compiled in the table 1 on working paper in overlead.
The files are mostly from the path of the fold called ODE, file paths for parameter estimated of simulations are:

# Empirical data of Patella species:
"Lapacom/Data/ToAnalyze/Madeira/BD_LIMPETS_MAD_1996-2018.xlsx"

 # Genetic data for Patella ordinaria and Patella aspera in madeira sites: 
"Lapacom/Data/GeneticDistance/Genetic_data.docx"

# Vector of Morality rates for life states estimations (d_{i}):
"Lapacom/Code/julia/ODE/estimate_mortality_rates.jl"

## Distance matrix estimation for estimate migration probability matrix (p_{kl}):
# .jl file
"Lapacom/Code/julia/ODE/distance_matrix.jl"
# .csv file
"Lapacom/Code/julia/ODE/distance_matrix.csv"

# Migration probability matrix (p_{kl}):
"Lapacom/Code/julia/ODE/migration_probabilities_among_sites.csv

# Parameters of initial states of the ststems simulation of CLC.
"Lapacom/results/001/src/load_params.jl"

# Simulation of Complex Life Cycle for two species with competition:
"Lapacom/Code/julia/ODE/model_dev.jl"

# Simulations of Simple Life Cycle of two species with competition:
"Lapacom/Code/julia/ODE/Numerical_aproximation/X_scenarios.jl"

Specific parameters used in Simulations for Patella ordinaria and Patella aspera for empirical data: 
# Egg deposition average 
oocytes_po = 385613 # Average: Patella ordinaria (nº of Eggs)
oocytes_pa = 77404  # Average: Patella aspera (nº of Eggs)
# Competition rate estimations
c_po = oocytes_po/(oocytes_pa + oocytes_po) # Patella ordinaria
c_pa = oocytes_pa/(oocytes_po + oocytes_pa) # Patella aspera

oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
reggs = oocytes / (365 * 0.42)       # Population growth rate    
re = reggs*0.998611*0.971057*0.4820525*0.00629     # adult poplation growth rate vector considering previous life stages.
Kt = 64000          # Carrying capacity (Theoretical aproximation)
H_ = [0.639,0.57]   # Exploitation rate (H) from empirical reference gived in the workshops
da_ = [0.55,0.59]    # Natural mortality rate for adults from empirical reference  collected in the workshops 
Sm = 56              # Maximum size for adults from empirical reference collected in the workshop
gammas = [0.32,0.36] # Adult growth rate from empirical reference collected in the workshop
Naj_po = reggs[1] # Adult individual abundancec for competition simulations in Patella aspera conditions(cij * Na_i *Na_j) 
Naj_pa = reggs[2] # Adult individual abundancec for competition simulations in Patella ordinaria conditions
k_ = 0.42 # Reproduction time per year of Patella species. Empirical data collected in the workshops from empirical reference.
t0_ = 365.25 # Average amount of days over a year.

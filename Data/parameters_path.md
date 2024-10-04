# List of file paths 
The list of all parameter used in the simulation for Patella species life cycles are compiled in the table 1 on working paper in overleaf. The file paths for parameter estimated used in simulations are:

## Empirical data of Patella species:
"Lapacom/Data/ToAnalyze/Madeira/BD_LIMPETS_MAD_1996-2018.xlsx"

## Genetic data for Patella ordinaria and Patella aspera in madeira sites: 
"Lapacom/Data/GeneticDistance/Genetic_data.docx"

## Vector of Natual Morality rates for life states estimations (d_{i}):
"Lapacom/Code/julia/ODE/estimate_mortality_rates.jl"

## Distance matrix estimation for estimate migration probability matrix (p_{kl}):
### .jl file
"Lapacom/Code/julia/ODE/distance_matrix.jl"
### .csv file
"Lapacom/Code/julia/ODE/distance_matrix.csv"

## Migration probability matrix (p_{kl}):
"Lapacom/Code/julia/ODE/migration_probabilities_among_sites.csv

# Parameters of initial states of the ststems simulation of CLC.
"Lapacom/results/001/src/load_params.jl"

## Simulation of Complex Life Cycle for two species with competition:
"Lapacom/Code/julia/ODE/model_dev.jl"

## Simulations of Simple Life Cycle of two species with competition:
"Lapacom/Code/julia/ODE/Numerical_aproximation/X_scenarios.jl"

# Specific parameters used in Simulations for Patella ordinaria and Patella aspera from empirical data and as theoretical aproximations: 

## Average numbe of oocytes per spawning season (Empitical estimation from dempirical data.
oocytes_po = 385613 # Average: Patella ordinaria (nº of Eggs)

oocytes_pa = 77404  # Average: Patella aspera (nº of Eggs)

oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera

reggs = oocytes / (365 * 0.42)       # Population growth rate    
g_ij = [reggs, 0.998611, 0.971057, 0.4820525, 0.00629] #Instant conversion population growth rate between life stages for Patella species.

re = reggs*0.998611*0.971057*0.4820525*0.00629     # adult poplation growth rate vector considering previous life stages.

Kt = 64000          # Carrying capacity (Theoretical aproximation)

H_{i} = [0.639,0.57]   # Exploitation rate (H) from empirical reference gived in the workshops
da_ = [0.55,0.59]    # Natural mortality rates for Patella species adults from empirical data collected in the workshops 

Sm = 56              # Maximum size for adults from empirical data collected in the workshop

gammas = [0.32,0.36] # Adult growth rate from empirical data collected in the workshop

## Competition rate estimations Overleaf eq: 12

c_po = oocytes_po/(oocytes_pa + oocytes_po) # Patella ordinaria

c_pa = oocytes_pa/(oocytes_po + oocytes_pa) # Patella aspera

Complementary information for SLC simulation (Theoretical aproximation)
Naj_po = reggs[1] # Adult individual abundancec for competition simulations in Patella aspera conditions (cpa * Na_pa * Naj_po) 

Naj_pa = reggs[2] # Adult individual abundancec for competition simulations in Patella ordinaria conditions (cpo * Na_po * Naj_pa)

k_ = 0.42 # Reproduction time per year of Patella species from Empirical data collected in the workshops from empirical data.

t0_ = 365.25 # Average amount of days over a year.

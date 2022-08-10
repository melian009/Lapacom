## Lapacom
Metacommunity dynamics of complex life-cycles in exploited seascapes

## Questions

* Predicting extinction probabilities in two metapopulations with and without interactions along a exploitation-connectivity gradient 

* Spatiotemporal synchrony promoted by humans?

* Do control sites exhibit higher asyncrhony? (Madeira)

* Richness peaks in the sites with intermediate connections? (Lower richness in the more isolated (connected) sites?) :: Richness here can be the variance of a trait

## Data

* Switchdrive: https://drive.switch.ch/index.php/s/CnVXHOa0YHL2iMb
* Distance matrix: Empirical distance matrix with flow 
        
## Code 

* Github: https://github.com/melian009/Lapacom

* Preliminary Octave version of Leslie matrix in space 
* Julia updated version

## Modeling

# Theory
* ODEs system to explore the extinction probability in a exploitation-connectivity gradient in two scenarios: One random and the second heterogenous both with a fixed number of source-sinks ratios. 
* Language: julia
* Package: DifferentialEquation.jl

# Theory-empirical
* Fit the following scenarios to the empirical data for the case of Madeira:
* Language: julia
* Package: Flux.jl

#===========================================================================================================
# Open position for 5-6 months master student to bridge ecological data and modeling in exploited ecosystems

A 5-6 months master student position is open for the project Lapacom, 
Metacommunity dynamics of complex life-cycles in exploited seascapes, 
funded by the Fundacao para a Ciencia e a Tecnologia, Portugal. Lapacom aims at
connecting empirical data for intertidal species to extinction dynamics in exploited metacommunities. 
The candidate should have a burning interest in computing programming, i.e., mostly julia but also R and other
open-source programming languages. Main task involves building codes around empirical data under different
exploitation scenarios to predict extinction probabilities in two Limpet species living in the intertidal 
habitats of Madeira island. You will work with a theoretical and computational team lead by 
Dr. Ali Vahdati and Dr. Carlos Melian. The overall vision is to understand the impact of 
overexploitation on metacommunities composed by complex life-cycles species. 
Please contact X and/ or Y if you need further information about the position. 
Deadline for applications is on X. 
#============================================================================================================


## Model 1: Infinite sites in homogeneous landscapes
* Model a large number of sites with no differentiation between North and South of the island, along connectivity-exploitation gradient and symmetry-asymmetry dispersal matrix. 
* Assumptions:
  * biotic trait: weak competition
  * abiotic trait: environmental data 
  * mating trait: 
  * dispersal trait:

## Model 2: Infinite sites in heterogeneous landscapes 
* Model a large number of sites with differentitation between North and South, exposition to human settlements and asymmetry in dispersal probabilities.
* Assumptions:
  * biotic trait: weak competition
  * abiotic trait: 
  * mating trait: 
  * dispersal trait:

## Model 3: Finite sites in homogeneous landscapes 
* Model the sampled number of sites with differentiation between norte and south using a random distance matrix.
* Assumptions:
  * biotic: strong competition
  * abiotic: 
  * mating:
  * dispersal:

## Model 4: Finite sites in heterogeneous landscapes
* Model the sampled number of sites with differentiation between norte and south using a distance-currents-driven matrix.
* Assumptions:
  * biotic: strong competition
  * abiotic: 
  * mating:
  * dispersal: 


## WorkingPaper

* Overleaf: https://www.overleaf.com/4457345492bddrnnxtdpqs
* Draft outlining the intro 
* Check format for references

## References

* Mendeley: https://www.mendeley.com/reference-manager/library/groups/private/3d5a2b3a-4020-30bc-a6fb-d5b9a2386c3a/all-references/

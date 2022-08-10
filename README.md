## Lapacom
Metacommunity dynamics of complex life-cycles in exploited seascapes

## Questions

* Does complex life-cycle predit exploitation response compared to simple life-cycle?

* How do complex life-cycle alter extinction probabilities in two metapopulations with and without interactions along a exploitation-connectivity gradient?

* Do control sites exhibit higher asyncrhony? (Madeira)

* Richness peaks in the sites with intermediate connections? (Lower richness in the more isolated (connected) sites?) :: Richness here can be the variance of a trait

## Data

* Switchdrive: https://drive.switch.ch/index.php/s/CnVXHOa0YHL2iMb
* Distance matrix: Empirical distance matrix with flow 
        
## Code 

* Github: https://github.com/melian009/Lapacom

* Preliminary Octave version of Leslie matrix in space 
* Julia updated version -- TODO: extension to explore "Dispersal capacity decay with life-cycle (DCLC)" as a way to distinguish SLC from CLC.                                                                                                                                                    

## Modeling

# Theory
* ODEs system to explore the extinction probability in a exploitation-connectivity gradient in two scenarios accounting for complex life-cycle: One spatially random and the second heterogenous both with a fixed number of source-sinks ratios. 
* Language: julia
* Package: DifferentialEquation.jl

# Theory-empirical
* Fit distance to human vs number of reproductives and extinction probability under the SLC and CLC scenarios for the case of Madeira:
* Language: julia
* Package: Flux.jl

'''
Modeling scenarios: ODE system with SLC (linear dispersal capacity decay with life cycle (DCLC)) and with CLC (non-linear DCLC). There are four scenarios for each life-cycle:
'''

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

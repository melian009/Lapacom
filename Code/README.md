Four main scenarios -------------------------------------------------------------------------

Model 1: Infinite sites in homogeneous landscapes
This scenario represents a large number of sites with no differentiation between North and South, different rates of human-driven disturbance or asymmetry in dispersal probabilities. Details of the biotic (weak competition), abiotic, mating and dispersal traits

Model 2: Infinite sites in heterogeneous landscapes 
A large number of sites with differentitation between North and South, exposition to human settlements and asymmetry in dispersal probabilities
Details of the biotic (weak competition), abiotic, mating and dispersal traits

Model 3: Finite sites in homogeneous landscapes 
Details of the biotic (strong competition), abiotic, mating and dispersal traits

Model 4: Finite sites in heterogeneous landscapes 
Details of the biotic (strong competition), abiotic, mating and dispersal traits
---------------------------------------------------------------------------------------------


Pseudocode core ------------------------------------------------------------------------------------------------------

1. Spatial structure: static
   1.1 Distance matrix containing size of each site, exposition (currents), and extraction
2. Agents: spatiotemporal dynamics 
   2.1 Life cycle
       Egg(A);Trocophore(B);Veliger(C);Juvenile(D), and Adult(E)
       --->>>Transition probabilities between stages driven by #post-fertilization hours  -->> Numbers in folder 'PPT_Outlines of development.pptx'
   2.2 Traits

       Mating 
          External fecundation --->> explicit males and females
          Random encounters
          Exponential decay #fertilized eggs per individual per site 

       Abiotic
          Exposition(North (less accessible) or South(more accessible))-extraction trade-offs 
          Currents to calculate between site distance

       Biotic
          Close season (Madeira)
          Extraction quantity:
                         Madeira:40mm or 15 Kg max
                         Canarias:45mm -- 10Kg professionals and 3kg non-professionals
          Distance to human settlements (near-far)
          Extractions independent of size
          Probability to die independent of size

       Dispersal
          Larger size more eggs (30-35mm first maturity size :: fertile)
          Global dispersal (A,B)
          Regional dispersal (B,C)
          Local dispersal (C)
          No dispersal (D,E)
          

3. Outputs

   Immigrant-Emigrants ratio ---->> sources and sinks 
   MPA (Marine protected areas) overlap sources and sinks?
   Local extinction probability each scenario :: 
   Local extinction is when number individuals 35mm equal to 0
   Plot x-axis distance to human-settlements vs y-axis %inds >35mm:
   Which is the scenario that best predict the empirical pattern?



4. ABC :: Determine priors y epsilon(posteriors) for the different parameters 


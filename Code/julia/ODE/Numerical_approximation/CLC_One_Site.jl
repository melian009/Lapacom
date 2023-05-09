## Packeges
using Pkg
Pkg.activate(".")
using LinearAlgebra
# using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim


#=
Formulation of 3 types of  life cicle for one site:

- Simple Life Cicle:
    Population dynamic model:

    dNa/dt = X * r * Ne * ((K - Na)/K) - da*Na - H * (1 - X) * Na
    
- Complex Life Cicle:

    Population dynamic models:

    - Two stages: Larvae/Egg & Adults.

    dNe/dt = X * r * Na - gEA * Na - de * Na
    dNa/dt = r * Na * ((K-Na)/K) - da * Na - (1-X) * H * Na

    - Five stages: Eggs, Veliger, Trocophore, Juvenile & Adults.

    dNe/dt = X * r * Ne - get * Ne - de * Ne
    dNt/dt = get * Ne - gtv * Nt - dt * Ne
    dNv/dt = gtv * Nt * ((K-Nv)/K) - gvj * Nv - dv * Nv
    dNj/dt = gvj * Nv * ((K-Nj)/K) - gja * Nj - dj * Nj
    dNa/dt = gja * Nj * ((K-Na)/K) - da * Na - H * (1 - X) * Na

    Metapopulation dynamic models 
    - Two stages: Larvae/Egg & Adults.

    dNa/dt = X * r * Na - gEA * Na - de * Na + sum(pkl*Nek) - Nel*sum(plk)
    dNa/dt = r * Na * ((K-Na)/K) - da * Na - (1-X) * H * Na

    - Five stages: Eggs, Veliger, Trocophore, Juvenile & Adults.

    dNe/dt = X * r * Ne - get * Ne - de * Ne + sum(pkl*Nek) - Nel*sum(plk)
    dNt/dt = get * Ne - gtv * Nt - dt * Ne + sum(pkl*Ntk) - Ntl*sum(plk)
    dNv/dt = gtv * Nt * ((K-Nv)/K) - gvj * Nv - dv * Nv + sum(pkl*Nvk) - Nvl*sum(plk)
    dNj/dt = gvj * Nv * ((K-Nj)/K) - gja * Nj - dj * Nj
    dNa/dt = gja * Nj * ((K-Na)/K) - da * Na - H * (1 - X) * Na

 
 
 
 
 
 
    Metacommunity dynamic model:

    dNe(i)/dt = X * r(i) * Na(i) - gEA * Na(i) - de(i) * Ne(i)
    dNa(i)/dt = gEA * Na(i) * ((K - (Na(i) + Na(i + 1))/K) - dA(i) * Na(i) - (1-X) * H(i) * Na(i)
    dSa(i)/dt = gamma * Sa(i) * (1 - (Sa(i)/(Smax-Smax * H(i) * X)))
    
    where i = nº of species == ["Patella ordinaria", "Patella aspera"]

Parameters and variables:
 - Ne = eggs abundance
 - Na = adults abuncance
 - Sa = adults size (Average sizes Before MPS+FULL = [46.44,44.45])
 - r = population growth rate [9.17,5.03]
 - K = carrying capacity (k =1e^4)
 - X = Reproductive period [1,0] 
 - (1-X) = Exploitation periosd [1,0]
 - H = Exploitation rate (H = [0.639,0.57])
 - gEA = instant conversion rate of life stages (EA = Eggs to Adults) (gEA = 0.006)
 - de = natural mortality rate or death rate for eggs (de = [de_po,de_pa]) # Note: this value needs to be defined.
 - da = natural mortality rate or death rate for adults (da = [0.55,0.59]) # Note: empirical estimated values
 - Smax = maximum adult size estimated (56.0mm) # Note: Empirical value from the sample analized
 - gamma = adult growth rate (gamma=[0.32,0.36] year^{-1})

Empirical estimated mortalities (Z,d,F) and exploitation rates (E):

 {Patella ordinaria} (Henriques, 2012)

 - Natural Mortality rate (d) was 0.55 year^{-1}
 - Fishing mortality rate (F) was 1.24 year^{−1} 
 - Total mortality rate (Z=d+F) was 1.79 year^{−1}; 
 - Exploitation rate (E=F/Z) was 0.693. 

 {Patella aspera} (Sousa, 2017)

 - Natural mortality (d) was 0.59 year^{-1};
 - Fishing mortality (F) was 0.79 year^{-1};
 - Total mortality (Z=d+F) was 1.38 year^{-1};
 - Exploitation rate (E=F/Z) was 0.57. 

 For numerical simulation use these rates (d and E).

Average sizes before and after marine protected area implementations

 Before (1996-2006): FULL ACCESS.
 Patella apera = 43.53mm
 Patella ordinaria = 46.26mm

 After  (2007-2017): MPA+FULL
 Patella aspera = 44.45mm
 Patella ordinaria = 46.44mm

 Only MPS (2007-2017)
 Patella aspera = 50.61mm
 Patella ordinaria = 49.25mm

 Only Full Access (2007-2017) 

 Patella aspera = 43.41mm
 Patella ordinaria = 45.72mm
=#



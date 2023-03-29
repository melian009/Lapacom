#Packages
using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Symbolics
using Plots
using Symbolics
import ForwardDiff.jacobian

#=
Formulation of the simple life cicle for one site:

    Metapopulation dynamic model:

    dNe/dt = X * r * Na - gEA * Na - de * Na
    dNa/dt = gEA * Ne * ((K-Ne)/K) - da * Na - (1-X) * H * Na
    dSa/dt = gamma * Sa * (1-(Sa/(Smax - Smax * H * X)))

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

# Metapopulation dynamic model 

function SLC_metapop!(du, u, p, t)
    Ne, Na, Sa = u
    i, r, K, rate, Exp, X, g, de,da, Smax, gamma = p
    du[1] = dNe = (X(t) * r[i] * Na) - (de[i] * Ne) - (g * Na)
    du[2] = dNa = (g * Ne * ((K - Na) / K)) - (da[i] * Na) - (Exp(t,rate,i) * Na)
    du[3] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 * Exp(t,rate,i))))
end

#Exploitation period/rate implamantation H(t,rate,i)= (1-X)*H(i)

function H(t, rate, i) #(1-X)*H(i)
    if (t % 365) / 365 < 0.42
      return rate[i]
    else
      return 0.0
    end
end

# Reproductive Cycleç

function rep(t)
    if (t % 365) / 365 >= 0.42
      return 1.0
    else
      return 0.0
    end
end

# Population Growth rate estimation (r=reggs):

oocytes_po = 385613 + 194902 - 194902  # Average: Patella ordinaria (nº of Eggs)
oocytes_pa = 703029 + 43496 - 43496  # Average: Patella aspera (nº of Eggs)
oocytes = [oocytes_pa,oocytes_po]    # Patella ordinaria, Patella aspera
reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.

re = reggs / 500     # Population growth rate
Kt = 1e4             # Carrying capacity
rates = [0.639,0.57] # Exploitation rate
gEA = 0.006          # Instant conversion between stages.
da_ =  [0.55,0.59]   # Natural mortality rate for adults
de_ = [0.975,0.977]  # Not estimated values. Need to be calculated by numerical aproximation
Sm = 56              # Maximum size for adults
gammas = [0.32,0.36] # Adult growth rate
i = [1,2]            # Species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)

P_sol_po = [i[1], re, Kt, rates, H, rep, gEA, de_,da_, Sm, gammas] # "Patella ordinaria" 
P_sol_pa = [i[2], re, Kt, rates, H, rep, gEA, de_,da_, Sm, gammas] # "Patella aspera" 

u0 = [1e3,1e3,33.4; 1e3,1e3,34.6]     # Initial populations abundance and size for the simulations

tspan = (0.0,365.0*2) # Temporal ranges for simulations: 2 years.

#Simulation for "Patella ordinaria"
prob_po = ODEProblem(SLC_metapop!, u0[1,], tspan, P_sol_po) 
sol_po = solve(prob_po, Tsit5())
# Simulation for "Patella aspera"
prob_pa = ODEProblem(SLC_metapop!, uo[2,], tspan, P_sol_pa)
sol_pa = solve(prob_pa,Tsit5())

Plots.plot(sol_1_po,vars=(0,1), label="Ne")
Plots.plot!(time,sol_1_po[2], label="Na")
Plots.title!("3rd definition: 1y")
Plots.xlabel!("t (days)")
Plots.ylabel!("N (Nº individuals)")





#Numerical aproximation to obtain natural mortality rates of eggs from empirical parameters.

@variables de i, r, K, rate, Exp, X, g,da, Smax, gamma, Ne, Na, Sa

dNe = (X(t) * r[i] * Na) - (de[i] * Ne) - (g * Na)
dNa = (g * Ne * ((K - Na) / K)) - (da[i] * Na) - (Exp(t,rate,i) * Na)
dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 * Exp(t,rate,i))))


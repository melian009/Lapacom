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
function SLC_metapop_before!(du, u, p, t)
  Ne, Na, Sa = u
   i, r, K, H, X, g, de,da, Smax, gamma = p
  du[1] = dNe = (X(t) * r[i] * Na * ((K - Na) / K)) - (de[i] * Ne) - (g * Ne)
  du[2] = dNa = (g * Ne * ((K - Ne) / K)) - (da[i] * Na) - (H[i] * Na)
  du[3] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 - H[i])))
end

function SLC_metapop_FULL!(du, u, p, t)
    Ne, Na, Sa = u
     i, r, K, H, X, g, de,da, Smax, gamma = p
    du[1] = dNe = (X(t) * r[i] * Na * ((K - Na) / K)) - (de[i] * Ne) - (g * Ne)
    du[2] = dNa = (g * Ne * ((K - Ne) / K)) - (da[i] * Na) - ((1 - X(t))* H[i] * Na)
    du[3] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 * (1 - X(t)) * H[i])))
end

function SLC_metapop_MPA!(du, u, p, t)
  Ne, Na, Sa = u
   i, r, K, X, g, de,da, Smax, gamma = p
  du[1] = dNe = (X(t) * r[i] * Na * ((K - Na) / K)) - (de[i] * Ne) - (g * Ne)
  du[2] = dNa = (g * Ne * ((K - Ne) / K)) - (da[i] * Na)
  du[3] = dSa = gamma[i] * Sa * (1 - Sa / (Smax))
end

#Exploitation period/rate implamantation  X=0

# Reproductive Cycle (X=1)

function rep(t)
  if (t % 365) / 365 >= 0.42
    return 1.0 # Reproductive Cycle
  else
    return 0.0 # Exploitation Cycle
  end
end


# Population Growth rate estimation (r=reggs):

oocytes_po = 385613 + 194902 - 194902  # Average: Patella ordinaria (nº of Eggs)
oocytes_pa = 703029 + 43496 - 43496  # Average: Patella aspera (nº of Eggs)
oocytes = [oocytes_pa,oocytes_po]    # Patella ordinaria, Patella aspera
reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.

re = reggs / 500     # Population growth rate
Kt = 64000           # Carrying capacity
rates = [0.639,0.57] # Exploitation rate (H)
rates2 = [0.02,0.01]
gEA = 0.006          # Instant conversion between stages.
da_ = [0.55,0.59]    # Natural mortality rate for adults
de_ = [0.001,0.003]  # Not estimated values. Need to be calculated by numerical aproximation
Sm = 56              # Maximum size for adults
gammas = [0.32,0.36] # Adult growth rate
i = [1,2]            # Species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)


# Before and after
# Full acces:     i,    r,  K,  H,     X,   g,   de,  da,  Smax, gamma = p
P_sol_po_full = [i[1], re, Kt, rates, rep, gEA, de_, da_, Sm, gammas] # "Patella ordinaria" 
P_sol_pa_full = [i[2], re, Kt, rates, rep, gEA, de_, da_, Sm, gammas] # "Patella aspera" 
# After
#MPA:            i,    r,  K,  X,   g,   de,  da, Smax, gamma = p
P_sol_po_mpa = [i[1], re, Kt, rep, gEA, de_, da_, Sm, gammas] # "Patella ordinaria" 
P_sol_pa_mpa = [i[2], re, Kt, rep, gEA, de_, da_, Sm, gammas] # "Patella aspera" 


# Initial populations abundance and mean size  for the simulations
# Before:
# Full access (1996-2006)
u0_po_before = [1e4, 1e4, 46.26]    # Patella ordinaria 
u0_pa_before = [1e4, 1e4, 43.53]    # Patella aspera

# After:
# FULL + MPA  (2007-2017)
u0_po = [1e4, 1e4, 33.4]    # Patella ordinaria 
u0_pa = [1e4, 1e4, 34.6]    # Patella aspera

#Full access (2007-2017)
u0_po_full = [1e4, 1e4, 43.41]    # Patella ordinaria 
u0_pa_full = [1e4, 1e4, 45.72]    # Patella aspera

#MPA (2007-2017
u0_po_mpa = [1e4, 1e4, 49.25] #Patella ordinaria 
u0_pa_mpa = [1e4, 1e4, 50.61] #Patella aspera


tspan = (1, 3000) # Temporal ranges for simulations: 2 years.

#Simulation for "Patella ordinaria"

#Before
prob_po_before = ODEProblem(SLC_metapop_FULL!, u0_po_before, tspan, P_sol_po_full) 
sol_po_before = solve(prob_po_before, Tsit5())

#After:
#Full Acces
prob_po_full = ODEProblem(SLC_metapop_FULL!, u0_po_full, tspan, P_sol_po_full) 
sol_po_full = solve(prob_po_full, Tsit5())

#MPA
prob_po_mpa = ODEProblem(SLC_metapop_MPA!, u0_po_mpa, tspan, P_sol_po_mpa) 
sol_po_mpa = solve(prob_po_mpa, Tsit5())


# Simulation for "Patella aspera"
#Before:
prob_pa_before = ODEProblem(SLC_metapop_before!, u0_pa_before, tspan, P_sol_pa_full)
sol_pa_before = solve(prob_pa_before,Tsit5())

# After:
# Full Acces
prob_pa_full = ODEProblem(SLC_metapop_before!, u0_pa_full, tspan, P_sol_pa_full)
sol_pa_full = solve(prob_pa_full,Tsit5())

# MPA
prob_pa_mpa = ODEProblem(SLC_metapop_MPA!, u0_pa_mpa, tspan, P_sol_pa_mpa)
sol_pa_mpa = solve(prob_pa_mpa,Tsit5())




#Plots
#Patella ordinaria
plot(sol_po_mpa, vars=(0,1), yscale=:log10,  label= "Ne: Before")
plot!(sol_po_mpa, vars=(0,1), yscale=:log10,  label= "Ne: Full")
plot!(sol_po_mpa, vars=(0,1), yscale=:log10,  label= "Ne: MPA")
plot!(sol_po_before, vars=(0,2), yscale=:log10, label= "Na: Before")
plot!(sol_po_full, vars=(0,2), yscale=:log10, label= "Na: Full access")
plot!(sol_po_mpa, vars=(0,2), yscale=:log10,  label= "Na: MPA")
title!("'Patella ordinaria'")
xlabel!("t (days)")
ylabel!("LOG10(N) (Nº individuals)")

#Patella ordinaria
plot(sol_pa_before, vars=(0,2), yscale=:log10, label= "Na: Before")
plot!(sol_pa_full, vars=(0,2), yscale=:log10, label= "Na: Full access")
plot!(sol_pa_mpa, vars=(0,2), yscale=:log10,  label= "Na: MPA")
title!("'Patella aspera'")
xlabel!("t (days)")
ylabel!("LOG10(N) (Nº individuals)")

plot(sol_pa_before, vars=(0,2), yscale=:log10, label= "Na: Before")
plot!(sol_pa_full, vars=(0,2), yscale=:log10, label= "Na: Full access")
plot!(sol_pa_mpa, vars=(0,2), yscale=:log10,  label= "Na: MPA")
title!("'Patella aspera'")
xlabel!("t (days)")
ylabel!("LOG10(N) (Nº individuals)")





# Numerical aproximation: jacobian matrix determination.

@variables Na, Ne , Sa, r, K, rate, Exp, X, g, de, da, Smax, gamma

# Symbolics.jacobian([f1(y1,y2), f2(y1,y2)],[y1, y2])

J = Symbolics.jacobian([(X * r * Na) - (de * Ne) - (g * Ne),
 (g * Ne * ((K - Ne) / K)) - (da * Na) - ((1 - X) * rate * Na), 
 gamma * Sa * (1 - Sa / (Smax * (1 - rate*(1-X))))],
 [Na])


# Cálculo del determinante
Det_J = det(J)

#Simplificación del determinante
M = Symbolics.simplify(Det_J)


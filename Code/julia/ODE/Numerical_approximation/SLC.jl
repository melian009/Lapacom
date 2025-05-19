#Packages
using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using GlobalSensitivity
#using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Symbolics
using GLMakie
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
  Na, Sa = u
   i, r, K, H, X, g, da, Smax, gamma = p
  #du[1] = dNe = (X(t) * r[i] * Na * ((K - Na) / K)) - (de[i] * Ne) - (g * Ne)
  du[2] = dNa = (g * Ne * ((K - Ne) / K)) - (da[i] * Na) - (H[i] * Na)
  du[3] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 - H[i])))
end

function SLC_metapop_FULL!(du, u, p, t)
    Na, Sa = u
     i, r, K, H, X, g, da, Smax, gamma = p
    #du[1] = dNe = (X(t) * r[i] * Na * ((K - Na) / K)) - (de[i] * Ne) - (g * Ne)
    du[2] = dNa = (g * Ne * ((K - Ne) / K)) - (da[i] * Na) - ((1 - X(t))* H[i] * Na)
    du[3] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 * (1 - X(t)) * H[i])))
end

function SLC_metapop_MPA!(du, u, p, t)
  Na, Sa = u
   i, r, K, X, g, de,da, Smax, gamma = p
  #du[1] = dNe = (X(t) * r[i] * Na * ((K - Na) / K)) - (de[i] * Ne) - (g * Ne)
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

oocytes_po = 385613                  # Average: Patella ordinaria (nº of Eggs)
oocytes_pa = 77404                  # Average: Patella aspera (nº of Eggs)
oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.

re = reggs / 500     # Population growth rate
Kt = 640000          # Carrying capacity
rates = [0.639,0.57] # Exploitation rate (H)
rates2 = [0.02,0.01]
gEA = 0.006          # Instant conversion between stages.
da_ = [0.55,0.59]    # Natural mortality rate for adults
Sm = 56              # Maximum size for adults
gammas = [0.32,0.36] # Adult growth rate
i = [1,2]            # Species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)

#de_ = [de_po,de_pe] # Not estimated values. Need to be calculated by numerical aproximation

#= Numerical aproximation: jacobian matrix determination for Full acces SLC metapopulation model

@variables Na, Ne , Sa, r, K, rate, Exp, X, g, de, da, Smax, gamma

# Symbolics.jacobian([f1(y1,y2), f2(y1,y2)],[y1, y2])

J = Symbolics.jacobian([(X * r * Na) - (de * Ne) - (g * Ne),
 (g * Ne * ((K - Ne) / K)) - (da * Na) - ((1 - X) * rate * Na), 
 gamma * Sa * (1 - (Sa / (Smax * (1 - rate*(1-X)))))],
 [Sa])
 expand(J)

# Cálculo del determinante
Det_J = det(J)
#Simplificación del determinante
M = Symbolics.simplify(Det_J)

# Depejamos Ne del Det_J sabiendo que para que se produzcan huevos X = 1.
Ne = 

Ne = K/2*(1-((de*da+g*de)/r*g))

#Despejamos de del determinante

de=((2*K*g)/(K*(da+g)))*((K/2)-Ne)
=#

# En base a la cantidad de huevos promedio que pone cada especie,
# definimos una capacidad de carga de huevos común para ambas metapoblaciones (KTotal).

KTotal=oocytes_pa+oocytes_po

# Para "Patella ordinaria".
# Calculamos la cantidad de huevos que mueren de una poblacion inicial promedio.

Ne_po=oocytes_po

dNe_po=((2*KTotal*gEA)/(KTotal*(da_[1]+gEA)))*((KTotal/2)-Ne_po)

# Calculamos la proporción de huevo muertos en relación al total de huevos iniciales que consideramos.
de_po =abs(dNe_po/Ne_po)

# Repetimos el mismo procedimiento para la otra especie: "Patella asera".

Ne_pa=oocytes_pa

dNe_pa=((2*KTotal*gEA)/(KTotal*(da_[2]+gEA)))*((KTotal/2)-Ne_pa)

de_pa =abs(dNe_pa/Ne_pa)

# Definimos un vector que engloba ambos ratios de mortalidad de huevos para el periodo reproductivo.

de_=(de_po,de_pa)



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
prob_po_before = ODEProblem(SLC_metapop_before!, u0_po_before, tspan, P_sol_po_full) 
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
prob_pa_full = ODEProblem(SLC_metapop_FULL!, u0_pa_full, tspan, P_sol_pa_full)
sol_pa_full = solve(prob_pa_full,Tsit5())

# MPA
prob_pa_mpa = ODEProblem(SLC_metapop_MPA!, u0_pa_mpa, tspan, P_sol_pa_mpa)
sol_pa_mpa = solve(prob_pa_mpa,Tsit5())




#Plots
#Patella ordinaria
lines(sol_po_before, vars=(0,1), yscale=:log10,  label= "Ne (1996-2006)")
lines!(sol_po_before, vars=(0,2), yscale=:log10,  label= "Na (1996-2006)")
lines!(sol_po_full, vars=(0,1), yscale=:log10,  label= "Ne (Full access)")
lines!(sol_po_full, vars=(0,2), yscale=:log10, label= "Na (Full access)")
lines!(sol_po_mpa, vars=(0,1), yscale=:log10, label= "Ne (MPA)")
lines!(sol_po_mpa, vars=(0,2), yscale=:log10,  label= "Na (MPA)")
title!("'Patella ordinaria'")
xlabel!("t (days)")
ylabel!("LOG10(N) (Nº individuals)")
savefig!("SLC_po_N_Before_Full_MPA_.png")

lines(sol_po_before, vars=(0,3), label= "Sa (1996-2006)")
lines!(sol_po_full, vars=(0,3), label= "Sa (Full access)")
lines!(sol_po_mpa, vars=(0,3),  label= "Sa (MPA)")
title!("'Patella ordinaria'")
xlabel!("t (days)")
ylabel!("Sa (mm)")
ylims!(55,57)
savefig!("SLC_po_S_Before_Full_MPA_.png")


#Patella ordinaria
lines(sol_pa_before, vars=(0,1), yscale=:log10,  label= "Ne (1996-2006)")
lines!(sol_pa_before, vars=(0,2), yscale=:log10,  label= "Na (1996-2006)")
lines!(sol_pa_full, vars=(0,1), yscale=:log10,  label= "Ne (Full access)")
lines!(sol_pa_full, vars=(0,2), yscale=:log10, label= "Na (Full access)")
lines!(sol_pa_mpa, vars=(0,1), yscale=:log10, label= "Ne (MPA)")
lines!(sol_pa_mpa, vars=(0,2), yscale=:log10,  label= "Na (MPA)")
title!("'Patella aspera'")
xlabel!("t (days)")
ylabel!("LOG10(N) (Nº individuals)")
savefig!("SLC_pa_N_Before_Full_MPA_.png")


lines(sol_pa_before, vars=(0,3), label= "Sa: Before")
lines!(sol_pa_full, vars=(0,3), label= "Sa: Full access")
lines!(sol_pa_mpa, vars=(0,3),  label= "Sa: MPA")
title!("'Patella aspera'")
xlabel!("t (days)")
ylabel!("Sa (mm)")
ylims!(55,56.5)
savefig!("SLC_pa_S_Before_Full_MPA_.png")





lines(sol_po_full, vars=(0,2), label= "Full access: 'Patella ordinaria")
lines!(sol_pa_full, vars=(0,2),  label= "Full access: 'Patella aspera")
ylims!(0,2000)
savefig!("SLC_Adults_Full_access.png")


lines(sol_po_mpa, vars=(0,2),  label= "MPA: 'Patella ordinaria")
ines!(sol_pa_mpa, vars=(0,2),  label= "MPA: 'Patella aspera'")
ylims!(0,2000)
savefig!("SLC_Adults_MPS.png")

lines(sol_po_mpa, vars=(0,2),  label= "MPA: 'Patella ordinaria")
lines!(sol_pa_mpa, vars=(0,2),  label= "MPA: 'Patella aspera'")
ylims!(0,2000)
savefig!("SLC_Adults_Full_access.png")
#Abuncances of spp A vs Abundances of spp B





#
#Exp_lim = 0.9999                 # Exploitation max limit 
#m = 0.0001                       # Interval of exploitation values
 
for n in 0:m:Exp_lim
function my_ode!(du, u, t, p)
    x1, x2, y1, y2 = u
    re1, re2, K, H1, H2, X, R1, R2, g1, g2, d1, d2, c12, c21, Smax = p
    
    #Parameter values species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)
    # Population Growth rate estimation (r=reggs):

   # oocytes_po = 385613                  # Average: Patella ordinaria (nº of Eggs)
   # oocytes_pa = 73029                   # Average: Patella aspera (nº of Eggs)
   # oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
   # reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.
   # re = reggs / 500     # Population growth rate: re[1,1] re[2,1]
    
   R1 = 1
   R2 = 1
   X =  0
   K = 640000          # Carrying capacity 
   H1 = n     #0.639
   H2 = 0.57 # Exploitation rate (H)
   g1 = 0.32
   g2 = 0.36 # Adult growth rate
    #rates2 = [0.02,0.01]
    #gEA = 0.006          # Instant conversion between stages.
   d1 = 0.55
   d2 = 0.59    # Natural mortality rate for adults
   Smax = 56              # Maximum size for adults
   re1 = 0.32
   re2 = 0.36 # Adult growth rate
   c12 = 0.05#competition term species 2 on 1
   c21 = 0.05#competition term species 1 on 2
    
    du[1] = (X * re1 * R1 * x1)*((K - x1)/K) - (d1 * x1) - H1*(1 - X)*x1 - c12*x1*x2
    du[2] = (X * re2 * R2 * x2)*((K - x2)/K) - (d2 * x2) - H2*(1 - X)*x2 - c21*x1*x2
    du[3] = (g1 * y1) * (1 - (y1/(Smax * (1 - H1 * (1-X)))))
    du[4] = (g2 * y2) * (1 - (y2/(Smax * (1 - H2 * (1-X)))))
end

u0 = [100.0, 100.0, 25.0, 25.0]  # Initial conditions
tspan = (0.0, 1000.0)  # Time span for the simulation (from t=0 to t=1000)
solver = Tsit5()
prob = ODEProblem(my_ode!, u0, tspan)
sol = solve(prob, solver)

plot!(sol, xlabel="Time", ylabel="State Variables", label=["x1" "x2" "y1" "y2"])

end

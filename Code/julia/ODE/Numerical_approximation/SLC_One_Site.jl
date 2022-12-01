## Packeges
using Pkg
Pkg.activate(".")
using LinearAlgebra
#using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
#using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Plots

### ------------------------------------------------------------------------------------------------
### 1. Model formulation for Simplex Life Cycle (SLC)  
### ------------------------------------------------------------------------------------------------
"""
Parameters:
  - r: intrinsic growth rate
  - g: instant conversion rate of eggs turning into adult
  - d: death rate
  - E: exploitation rate
  - K: carrying capacity
  - X_rc: Reproductive Cycle
"""
#= 
In the first definition of the equations, the abundances equatons arent linked with the size.
=#

# The model definition is:

function single_site!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = p
  du[1] = dNⱼ = (r * Nₐ * ((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
end

#=
In the second definition, the three equations are defined in such a way that
the number of eggs (N_e) depends on the size of the adults (S_a) and the adult abundance (N_a)
which depends on the exploitation rate (Et).
This means that the three ODE are linked between each other.
=#

# The ODE definition is:

function SLC_single_site_S!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = p
  du[1] = dNⱼ = (r * Nₐ * (Sₐ/sizeₘₐₓ)*((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
end

#= 
In this last definition, the equiations consider the reproductive/exploitation cycle (X_rc) that change during the year.
=#

# The ODE model definition is:

function SLC_single_site_S_X!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ, X_rc = p
  du[1] = dNⱼ = (X_rc * r * Nₐ * (Sₐ/sizeₘₐₓ)*((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - ((1 - X_rc) * E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ * (1 - E(t) * (1-X_rc))))
end

# The reproductive cycle is:

function X_rc(t)       
  if (t % 365) / 365 >= 0.42
    return 1.0        # equals to zero when fishing
  else
    return 0.0        # equals to one when adults reproduce
  end
end

#The time varing exploitation is:

Et(t) = (sin(t)^2) / 2  

# To adjust the sin function horizonally (stretch/shrink the wave length), multiply t by a factor c.

# and you may also use a function with if statements:

function Et(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return 0.5
  end
end

### -----------------------------------------------------------------------------------
### 2. Graphic representation of SLC in a single site for the three model formulations.
### -----------------------------------------------------------------------------------

# Fist definition: Abundances ODE and Size ODE are not linked and no reproductive cicle consideration:

# With function SLC_single_site!

# Parameters and initial conditions for the simulation 

p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0] # r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ 
u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ_0,Nₐ_0, Sₐ_0
tspan_1 = (0.0, 50.0)
tspan_2 = (0.0, 365.0)
prob_1 = ODEProblem(single_site!, u0_1, tspan_1, p_1)
sol_1 = solve(prob_1, Tsit5())

# Abundances plot in time 
plot(sol_1,vars=(0, 1), label="Nⱼ")
plot!(sol_1,vars=(0, 2), label="Nₐ")
title!("1st definition")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

# Adult size plot by time
plot(sol_1,vars=(0,3), label="Sₐ") 
title!("1st definition")
xlabel!("t (days)")
ylabel!("S (mm)")


"""
Second definition: Abundances ODE and Size ODE are linked, but no reproductive cicle consideration:
- with function SLC_single_site_S!
"""
#Parameters and initial conditions for the simulation 

p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0] # r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ 
u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ_0,Nₐ_0, Sₐ_0
tspan_1 = (0.0, 365.0)
prob_2 = ODEProblem(SLC_single_site_S!, u0_1, tspan_1, p_1) 
sol_2 = solve(prob_2, Tsit5())

#Abundances plot in time 
plot(sol_2,vars=(0, 1), label="Nⱼ")
plot!(sol_2,vars=(0, 2), label="Nₐ")
title!("Abundances for SLS in a sigle site: 1st definition")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

#Adult size plot by time
plot(sol_2,vars=(0,3), label="Sₐ") 
title!("Adult Size for SLS in a sigle site: 1st definition")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")













### ------------------------------------------------------------------------------
### 2. Numerical aproximation of SLC one site and no reproductive cicle inclusion
### ------------------------------------------------------------------------------

Exp_lim = 0.9999                 # Exploitation max limit 
m=0.0559                          # Interval of exploitation values 
Expl= 0:m:Exp_lim                # Expoitation values for plotting
N_et = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
N_at = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
S_at = zeros(Float64,size(Expl)) # Void vector to array the size of adults for diferent exploitation values
c = 0                            # C is the position of the vector N_at and N_et

for n = 0:m:Exp_lim
  #=
  function Et(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return n                    #For an exploitation value equal to 1, the mathematical result is erratic because the size equation will present a denominator division equal to 0
  end
 end
 =#
 Et(t) = (sin(t)^2) / 2  # time varing exploitation
 p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0]
 u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ_0,Nₐ_0, Sₐ_0
 tspan_1 = (0.0, 365.0) #time span range: One Year
 prob_1 = ODEProblem(single_site_SLC!, u0_1, tspan_1, p_1)
 sol_1 = solve(prob_1, Tsit5())
 c=c+1

 N_et[c,] = sol_1[1,end]
 N_at[c,] = sol_1[2,end]
 S_at[c,] = sol_1[3,end]

end


plot(Expl,N_et,label="Nⱼ")
plot!(Expl,N_at,label="Nₐ")
xlims!(0.0,1)
xlabel!("E")
ylabel!("N (nº individuals)")

plot(Expl,S_at)
xlabel!("E")
ylabel!("Sₐ (nº individuals)")




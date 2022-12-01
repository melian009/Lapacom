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
  - X: Reproductive Cycle
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
  r, g, dⱼ, dₐ, Exp, K, size_growth_rate, sizeₘₐₓ, X, Exp_rate = p
  du[1] = dNⱼ = (X(t) * r * Nₐ * (Sₐ/sizeₘₐₓ)*((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - ((1 - X(t)) * Exp(t,Exp_rate) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ * (1 - Exp(t,Exp_rate))))
end

# The reproductive cycle is:

function X_rc(t)
  if (t % 365) / 365 >= 0.42
    return 1.0
  else
    return 0.01
    # return 0.0
  end
end
 
# The time varing exploitation for simple model without consithering the reproductive cyceis:
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

#The time varing exploitation considerig the exploitation cicle is:

function Exp(t, Exp_rate)
  if (t % 365) / 365 < 0.42
    return Exp_rate
  else
    # return 0.0
    return 0.01
  end
end


### -----------------------------------------------------------------------------------
### 2. Graphic representation of SLC in a single site for the three model formulations.
### -----------------------------------------------------------------------------------

"""
Fist definition: Abundances ODE and Size ODE are not linked and no reproductive cicle consideration:
- With function SLC_single_site!
"""

# Parameters and initial conditions for the simulation 
p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0] # r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ 
u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ,Nₐ, Sₐ

# Calculate the numerical ODE solution for a 50 days time range
tspan_1 = (0.0, 50.0) 
prob_1 = ODEProblem(single_site!, u0_1, tspan_1, p_1)
sol_1 = solve(prob_1, Tsit5())

# Abundances plot in time 
plot(sol_1,vars=(0, 1), label="Nⱼ")
plot!(sol_1,vars=(0, 2), label="Nₐ")
title!("1st definition: 50 days")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

# Adult size plot by time
plot(sol_1,vars=(0,3), label="Sₐ") 
title!("1st definition: 50 days")
xlabel!("t (days)")
ylabel!("S (mm)")


# Calculate the numerical ODE solution for a 760 days (two years) time range
tspan_2 = (0.0, 365.0*2) 
prob_1_1 = ODEProblem(single_site!, u0_1, tspan_2, p_1)
sol_1_1 = solve(prob_1_1, Tsit5())


# Abundances plot in time 
plot(sol_1,vars=(0, 1), label="Nⱼ")
plot!(sol_1,vars=(0, 2), label="Nₐ")
title!("1st definition: 50 days")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

# Adult size plot by time
plot(sol_1,vars=(0,3), label="Sₐ") 
title!("1st definition: two years")
xlabel!("t (days)")
ylabel!("S (mm)")


"""
Second definition: Abundances ODE and Size ODE are linked, but no reproductive cicle consideration:
- With function SLC_single_site_S!
"""

#Parameters and initial conditions for the simulation are the same than the first definition

# Calculate the numerical ODE solution for a 365 days (a year) time range
prob_2 = ODEProblem(SLC_single_site_S!, u0_1, tspan_2, p_1) 
sol_2 = solve(prob_2, Tsit5())

#Abundances plot in time 
plot(sol_2,vars=(0, 1), label="Nⱼº")
plot!(sol_2,vars=(0, 2), label="Nₐ")
title!("2nd definition: two years")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

#Adult size plot by time
plot(sol_1_1,vars=(0,3), label="Sₐ:1º")
plot!(sol_2,vars=(0,3), label="Sₐ:2º")
title!("2nd definition")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")


"""
The third definition: Abundances ODE and Size ODE are linked and have the reproductive cicle consideration:
- With function SLC_single_site_S_X!
"""

#Parameters and initial conditions for the third simulation
p_2 = [0.6, 0.06, 0.05, 0.08, Exp, 1e4, 0.2, 40.0, X_rc, 0] # r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ, X, Exp_rate 
u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ,Nₐ, Sₐ


# Calculate the numerical ODE solution for a 365 days (a year) time range
tspan_3 =(0,365*2)
prob_3 = ODEProblem(SLC_single_site_S_X!, u0_1, tspan_3, p_2) 
sol_3 = solve(prob_3, Tsit5())

#Abundances plot in time 
plot(sol_3,vars=(0, 1), label="Nⱼ:3º")
plot!(sol_3,vars=(0, 2), label="Nₐ:3º")
title!("Abundances for SLS in a sigle site")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

#Adult size plot by time
plot(sol_1_1,vars=(0,3), label="Sₐ:1º")
plot!(sol_2,vars=(0,3), label="Sₐ:2º")
plot!(sol_3,vars=(0,3), label="Sₐ:3º")
title!("Adult Size for SLS in a sigle site")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")


### ------------------------------------------------------------------------------
### 2. Numerical aproximation of SLC one site and no reproductive cicle inclusion
### ------------------------------------------------------------------------------

""" ------------------------------------------------------------------------------
-                                    1st definition                              -
---------------------------------------------------------------------------------- """

Exp_lim = 0.9999                 # Exploitation max limit 
m=0.0559                          # Interval of exploitation values 
Expl= 0:m:Exp_lim                # Expoitation values for plotting
tspan = (0.0,365*2)              # Time range 
u0 = [1e3,1e3,40]                # Initial conditions of N_e, N_a, S_a

N_et_1 = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
N_at_1 = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
S_at_1 = zeros(Float64,size(Expl)) # Void vector to array the size of adults for diferent exploitation values
c = 0                              # C is the position of the vector N_et, N_at and S_at

for n = 0:m:Exp_lim
  function Et(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return n                    #For an exploitation value equal to 1, the mathematical result is erratic because the size equation will present a denominator division equal to 0
  end
 end
 p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0]
 prob_1 = ODEProblem(single_site!, u0, tspan, p_1)
 sol_1 = solve(prob_1, Tsit5())
 c=c+1

 N_et_1[c,] = sol_1[1,end]
 N_at_1[c,] = sol_1[2,end]
 S_at_1[c,] = sol_1[3,end]

end


plot(Expl,N_et_1,label="Nⱼ")
plot!(Expl,N_at_1,label="Nₐ")
xlims!(0.0,1)
xlabel!("E")
ylabel!("N (nº individuals)")

plot(Expl,S_at_1)
xlabel!("E")
ylabel!("Sₐ (nº individuals)")

""" ------------------------------------------------------------------------------
-                                    2nd definition                              -
---------------------------------------------------------------------------------- """
Exp_lim = 0.99999                # Exploitation max limit 
m=0.00001                        # Interval of exploitation values 
Expl= 0:m:Exp_lim                # Expoitation values for plotting
tspan = (0.0,365*2)              # Time range 
u0 = [1e3,1e3,40]                # Initial conditions of N_e, N_a, S_a


N_et_2 = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
N_at_2 = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
S_at_2 = zeros(Float64,size(Expl)) # Void vector to array the size of adults for diferent exploitation values
c = 0                              # C is the position of the vector N_et, N_at and S_at

for n = 0:m:Exp_lim
  function Et(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return n                    #For an exploitation value equal to 1, the mathematical result is erratic because the size equation will present a denominator division equal to 0
  end
 end
 p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0]
 prob_2 = ODEProblem(SLC_single_site_S!, u0, tspan, p_1)
 sol_2 = solve(prob_2, Tsit5())
 c=c+1

 N_et_2[c,] = sol_2[1,end]
 N_at_2[c,] = sol_2[2,end]
 S_at_2[c,] = sol_2[3,end]

end


plot(Expl,N_et_2,label="Nⱼ")
plot!(Expl,N_at_2,label="Nₐ")
xlims!(0.0,1)
xlabel!("E")
ylabel!("N (nº individuals)")


""" ------------------------------------------------------------------------------
-                                    3rd definition                              -
---------------------------------------------------------------------------------- """

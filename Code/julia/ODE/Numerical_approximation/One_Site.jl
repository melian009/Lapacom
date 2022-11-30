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

### ----------------------------------------------------------------
### 1. Simplex Life Cycle (SLC) in a single site
### ----------------------------------------------------------------
"""
r: intrinsic growth rate
g: instant conversion rate of eggs turning into adult
d: death rate
E: exploitation rate
K: carrying capacity
"""
#= In the first definition of the equations, the abundances equatons arent linked with the size one.

function single_site!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = p
  du[1] = dNⱼ = (r * Nₐ * ((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
end

In the new definition, the three equations are defined in such a way that the number of eggs (N_e)
depends on the size of the adults (S_a) and and the adult abundance which depends on the exploitation rate (Et).
This means that the three ODE are linked between each other.
=#
function single_site_SLC!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = p
  du[1] = dNⱼ = (r * Nₐ * (Sₐ/sizeₘₐₓ)*((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
end


Et(t) = (sin(t)^2) / 2  # time varying exploitation
# To adjust the sin function horizonally (stretch/shrink the wave length), multiply t by a factor c.
# You may also use a function with if statements.
function Et(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return 0.5
  end
end



p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0]
u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ_0,Nₐ_0, Sₐ_0
tspan_1 = (0.0, 365.0)
prob_1 = ODEProblem(single_site_SLC!, u0_1, tspan_1, p_1)
sol_1 = solve(prob_1, Tsit5())

plot(sol_1,vars=(0, 1), label="Nⱼ")
plot!(sol_1,vars=(0, 2), label="Nₐ")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")
plot(sol_1,vars=(0,3), label="Sₐ") #Size decrease in time from 40 to aproximate 27

### ----------------------------------------------------------------
### 2. Numerical aproximation of SLC
### ---------------------------------------------------------------- 

Exp_lim = 0.9999                 # Exploitation max limit 
m=0.0559                          # Interval of exploitation values 
Expl= 0:m:Exp_lim                # Expoitation values for plotting
N_et = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
N_at = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
S_at = zeros(Float64,size(Expl)) # Void vector to array the size of adults for diferent exploitation values
c = 0                            # C is the position of the vector N_at and N_et

for n = 0:m:Exp_lim
 function Et(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return n                    #For an exploitation value equal to 1, the mathematical result is erratic because the size equation will present a denominator division equal to 0
  end
 end
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
xlims!(0.0,0.6)
xlabel!("E")
ylabel!("N (nº individuals")

plot(Expl,S_at)
xlabel!("E")
ylabel!("Sₐ (nº individuals")



### ----------------------------------------------------------------
### Complex Life Cycle (CLC) in a single site
### ---------------------------------------------------------------- 

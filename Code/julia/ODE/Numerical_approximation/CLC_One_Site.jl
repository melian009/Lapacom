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

### ------------------------------------------------------------------------------------------------
### Model formulation for Complex Life Cycle (CLC) for a Single Site (One Site)
### ------------------------------------------------------------------------------------------------

"""
Variables:
  - N_e: Egg Abundance
  - N_t: Trochophore Abundance
  - N_v: Veliger Abundance
  - N_j: Juvenile Abundance
  - N_a: Adult Abundance
  - Sₐ: Adult Size
Parameters:
  - r: intrinsic growth rate
  - g: instant conversion rate of eggs turning into adult
  - d: death rate
  - E: exploitation rate
  - K: carrying capacity
  - X: Reproductive Cycle
"""

## Starting the model

#=    Simple life cycle: equations for one Site

    function single_site_S_X!(du, u, p, t)
      Nⱼ, Nₐ, Sₐ = u
      r, g, dⱼ, dₐ, Exp, K, size_growth_rate, size_max, X, rate = p
      du[1] = dNⱼ = (X(t) * r * Nₐ * (Sₐ/sizeₘₐₓ)*((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
      du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (Exp(t,rate) * Nₐ)
      du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ * (1 - Exp(t,rate))))
    end
=#

# Complex Life Cycle: equations for one Site

function CLC_OS!(du, u, p, t)
  N_e, N_t, N_v, N_j, N_a, S_a = u
  r, d, K, Exp, size_growth_rate, size_max, X, rate = p
  # dN_e/dt    = [X    * r    * N_a * (S/S_max)      * ((K - N_A)/k)] - (d_e * N_e)  - (g_et * N_e)
  du[1] = dN_e = (X(t) * r[1] * N_a * (S_a/size_max) * ((K * N_a)/K)) - (d[1] * N_e) - (r[2] * N_e)

  # dN_t/dt    = (g_et * N_e) - (d_t * N_t) - (g_tv * N_t)
  du[2] = dN_t = (r[2]* N_e) - (d[2] * N_t) - (r[3] * N_t)
  
  # dN_v/dt    = (g_tv * N_t) - (d_v * N_v) - (g_vj * N_v)
  du[3] = dN_v = (r[3]* N_t) - (d[3] * N_v) - (r[4] * N_v)

  #dN_j/dt     = (g_vj * N_v) - (d_j * N_j) - (g_ja * N_j)
  du[4] = dN_j = (r[4] * N_v) - (d[4] * N_t) - (r[5] * N_j)
  
  #dN_a/dt     = (g_ja * N_j) - (d_a * N_a) - [(1 - X) * E * N_a] 
  du[5] = dN_a = (r[5] * N_j) - (d[5] * N_a) - (Exp(t,rate) * N_a)

  #dS_a/dt     = γ * S_a * [1 - (S_a / (Smax * (1 - X)))]
  du[6] = dS_a = size_growth_rate * S_a * [1- S_a/(size_max * (1-Exp(t,rate)))]
end

## Model parameters definition:
# Exploitation: equals to zero when organisms reproduce

function Exp(t, rate)
  if (t % 365) / 365 < 0.42
    return rate
  else
    return 0.00
    end
end

# Captures the reproductive cycle are equals to zero when fishing and to one when organisms reproduce

function X(t)
  if (t % 365) / 365 >= 0.42
    return 1.0
  else
    return 0.00
  end
end


# Conversion rates between stages

# Average oocytes per year per adult

avg_oocytes = mean([92098, 804183]) #804.000 huevos para una talla 58mm 

# conversion rate of adults to eggs.

reggs = avg_oocytes / (365 * 0.42) 

reggs = reggs / 10000 

r=reggs #Eggs population intrinsic growth rate

# because the rate is too high to be handled
# = [r: 385.613eggs, g_et: 0.7d, g_tv: 1.3d, g_vj: 4d + 3d (7d), g_ja:730]
r = [reggs,          0.998611,   0.971057,   0.683772,           0.00629]

# natural death rates per life stage.
# = [d_e,   d_t,   t_v,   d_j,   d_a]
d = [0.001, 0.001, 0.001, 0.001, 0.000322]

size_growth_rate = 0.32/365 # γ 

exploitation_rates = rand(0.001:0.001:0.003)  # E: use empirical values

size_max = 56.0
K = 64_000  # for 6.4 km2 per site.

p_general =[r, d, K, Exp, size_growth_rate, size_max, X, exploitation_rates]

# Initial abundances of each life state and size.
u0 = [20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 40.0]
u0_ = convert(Float64, u0)

tspan_0 = (0,365*2) # Simulation for 2 years time range
tspan_0_ = convert(Tuple{Float64,Float64}, tspan_0)


prob_CLC_1= ODEProblem(CLC_OS!, u0, tspan_0_, p_general) 
sol_CLC_1 = solve(prob0_CLC_1, Tsit5())


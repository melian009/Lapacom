## Packeges
using Pkg
Pkg.activate(".")
using LinearAlgebra
#using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Plots
using Plots.PlotMeasures

### ------------------------------------------------------------------------------------------------
### 1. Model formulation for Complex Life Cycle (CLC) for a Single Site (One Site)
### ------------------------------------------------------------------------------------------------

"""
Variables:
  - Nⱼ: Egg Abundances
  - Nₐ: Adult Abundance
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
#---------------------------
#Nº of life stages

nlifestages = 5

## initial state of the system
# Initial population sizes at each life stage per site
# There are 100 individuals per m2. We assume 2 km2 -> 20k individuals.
u0_general = [1_800.0 for j in 1:nlifestages]

## defining model parameters
function exploit(t, rate)
  if (t % 365) / 365 < 0.42
    return rate
  else
    # return 0.0
    return 0.01
  end
end

"""
captures the reproductive cycle and equals to zero when fishing and to one when organisms reproduce
"""
function reproductive_cycle(t)
  if (t % 365) / 365 >= 0.42
    return 1.0
  else
    return 0.01
    # return 0.0
  end
end

# conversion rates between stages

# average oocytes per year per adult
avg_oocytes = mean([92098, 804183])

# conversion rate of adults to eggs.
reggs = avg_oocytes / (365 * 0.42) 
reggs = reggs / 500 
r=reggs
# because the rate is too high to be handled
 #r = [reggs, 0.998611, 0.971057, 0.683772, 0.00629]

# natural death rates per life stage.
d = [0.001, 0.001, 0.001, 0.001, 0.000322]

size_growth_rate = 0.32/365

exploitation_rates = rand(0.001:0.001:0.003)  # TODO: use empirical values

size_max = 56.0

K = 64_000  # for 6.4 km2 per site.

# α = [0.1, 0.1, 0.1]  # dispersion factor for Egg, Trochophore, and Veliger
# Since we do not have any info about site size, dispersion is only a function of dispersion factor and distance.
p_general =[r, d, size_growth_rate, exploitation_rates, size_max, K]


function Single_site_CLC!(du, u, p, t)
  # Change these parameters if you change the model
  #nsites = 4
  nlifestages = 5
  # site_indices = [0, 6, 12, 18, 24, 30, 36, 42]#, 48, 54]

  counter = 0
  #for site in 1:nsites
    # site_index = (site * (nlifestages + 1)) - (nlifestages + 1)
    # site_index = site_indices[site]

    # Stages 1 Egg
    counter += 1
    stage = 1
    prev_stage = 5
    du[counter] = (reproductive_cycle(t) * p[1][stage] * u[prev_stage] * ((p[6] - u[prev_stage]) / p[6])) -
                  (p[1][stage+1] * u[stage]) -
                  (p[2][stage] * u[stage])
    # Stage 2 Trochophore
    counter += 1
    stage = 2
    prev_stage = 1
    du[counter] = (p[1][stage] * u[prev_stage]) -
                  (p[1][stage+1] * u[stage]) -
                  (p[2][stage] * u[stage])
    #stage 3 Veliger
    counter += 1
    stage = 3
    prev_stage = 2
    du[counter] = (p[1][stage] * u[prev_stage]) -
              (p[1][stage+1] * u[stage]) -
              (p[2][stage] * u[stage])
    # stage 4 Juvenile
    counter += 1  
    stage = 4
    prev_stage = 3
    du[counter] = (p[1][stage] * u[prev_stage]) -
                (p[1][stage+1] * u[stage]) - 
                (p[2][stage] * u[stage]) 

    # stage 5 adult
    counter += 1
    stage = 5
    prev_stage = 4
    du[counter] = (p[1][stage] * u[site, prev_stage]) - 
                (exploit(t, p[4]) * u[ stage]) -
                (p[2][stage] * u[ stage])
                
    # adult sizes
    # adult sizes. dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
    counter += 1
    du[counter] = p[3] * u[nlifestages+1] * (1 - u[nlifestages+1] / (p[5] - (p[5] * exploit(t, p[4]))))
end

tspan_general = (100.0, 200.0)
prob_general = ODEProblem(Single_site_CLC!, u0_general, tspan_general, p_general)
# sol_general = solve(prob_general, Rosenbrock23());
# sol_general = solve(prob_general, alg_hints=[:stiff]);
sol_general = solve(prob_general, ROS34PW2(), dt=0.0000001, adaptive=false);

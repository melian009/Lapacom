using Pkg
Pkg.activate(".")
using LinearAlgebra
using OrdinaryDiffEq
# using DifferentialEquations
using GlobalSensitivity
using Plots
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim

## [x]: exploitation stops during certain months of the year. Implement time varying E.
## TODO: From topography, set distance between two sites to infinite, if there is land/another site between them
## [x]: Empirically estimate the rate at which individuals from one stage turn into another one:
## Conversion rates for the given number of days is as follows:
# 730: 0.00629
# 3: 0.784557
# 0.7: 0.998611
# 1.3: 0.971057
# 4: 0.683772
# In the ODEs, the unit of left side is individuals per day. This means that each time step is one day. The rate of the right side are calculated such that after the given time between life stages, 99% of the individuals from one stage are converted to the next stage. To that end, I solve the followin equation: 0.01 = r^t, where I replace t with the time between stages.
## [x]: Create another version of model where there are five life stages instead of two. In this system, the probability to migrate decreases exponentially. $1/d \times e^{-\gamma}$.
## [x]: Check how to add all equations in a loop instead of writing by hand.
# [x]: Modify the equations according to the following facts:
# 1. adults live 7 to 10 years, depending on their size. 
#   Age in days ranges between 2555 to 3650.
#   Actual age is (size * max age) / max size. If actual age < min age, actual age = min age. This is difficult to include in the ODE. We can ignore it and use an average age: 3102 days.
#   death probability of adults da =  1 / age. 1/3102 = 0.000322.
# 2. adults spawn 92,098 to 804,183 oocytes per year during the spawning season.
#   Appearance rate of eggs: #oocytes/365 in a sin function (to account for spawning season). Note that adults do not turn into eggs (they live many steps.)
## [x]: estimate `size_growth_rate` by knowing juvenile average size and adult average size and dividing their difference by 730 days. Growth rate is 0.32 per year, which is 0.32/365 per day.
## [ ] TODO: Estimate exploitation rates given average sizes and a single pool model.
#=
Example of adding if statements in the ODE system

function fun(du,u,p,t)
  if u[1] <30
      du[1] = (0.04*u[1] + 5)*u[1] + 150 - u[2] - p[5]
      du[2] = p[1]*(p[2]*u[1]-u[2])
  else
      u[1] = p[3]
      u[2] = u[2] + p[4]
  end
end
=#


#=
Five simplifying assumptions:

1. Transition times. Some short transition times can be aggregated.
2. The number of oocytes. We use the mean number of oocytes laid per adult per year.
3. Dispersal potential of trochophore and veliger. Eggs disperse passively, T and V can disperse actively.
4. Juveniles homing. Eggs have a global pool. Once they are released, they travel with the currents to the ocean. In the ocean they develop and actively migrate back to different sites. This is like the "lottery model".
5. Exploitation. Five months of exploitation, seven months rest. There is no reproduction during the exploitation time. NOTE that adults continueously keep growing.

=#

### ----------------------------------------------------------------
### 1. Single site
### ----------------------------------------------------------------

"""
r: intrinsic growth rate
g: rate of juveniles turning into adult
d: death rate
E: exploitation rate
K: carrying capacity
"""
function single_site!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = p
  du[1] = dNⱼ = (r * Nₐ * ((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
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
u0_1 = [1e3, 1e3, 40.0]
tspan_1 = (0.0, 50.0)
prob_1 = ODEProblem(single_site!, u0_1, tspan_1, p_1)
sol_1 = solve(prob_1, Tsit5())

plot(sol_1, vars=(0, 1), label="Nⱼ")
plot!(sol_1, vars=(0, 2), label="Nₐ")

### ----------------------------------------------------------------
### 2. Two sites
### ----------------------------------------------------------------

function two_sites!(du, u, p, t)
  Nⱼ₁, Nₐ₁, S₁,
  Nⱼ₂, Nₐ₂, S₂ = u

  r, g, dⱼ, dₐ, size_growth_rate,
  mₒᵤₜ₁, mₒᵤₜ₂,
  E₁, E₂,
  K₁, K₂,
  sizeₘₐₓ₁, sizeₘₐₓ₂ = p

  du[1] = dNⱼ₁ = (r * Nₐ₁ * ((K₁ - Nₐ₁) / K₁)) + (mₒᵤₜ₂ * Nⱼ₂) - (mₒᵤₜ₁ * Nⱼ₁) - (dⱼ * Nⱼ₁) - (g * Nⱼ₁)
  du[2] = dNₐ₁ = (g * Nⱼ₁) - (dₐ * Nₐ₁) - (E₁ * Nₐ₁)
  # track the influence of migration on body size. If there is immigration from a site where individuals are genetically small, then it should affect the distribution of body size in this site too. 
  # S₁ need to be replaced by the mean body size of this population which is a function of body sizes of all source populations and their contribution to this sink population.
  # frac = ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁)
  # du[3] = dS₁ = size_growth_rate * ((S₁ * (1 - ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁))) + (S₂ * ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁))) * (1 - ((S₁ * (1 - ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁))) + (S₂ * ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁))) / (sizeₘₐₓ₁ - (sizeₘₐₓ₁ * E₁)))
  du[3] = dS₁ = size_growth_rate * S₁ * (1 - S₁ / (sizeₘₐₓ₁ - (sizeₘₐₓ₁ * E₁)))  # logistic growth: rN(1- N/K)

  du[4] = dNⱼ₂ = (r * Nₐ₂ * ((K₂ - Nₐ₂) / K₂)) + (mₒᵤₜ₁ * Nⱼ₁) - (mₒᵤₜ₂ * Nⱼ₂) - (dⱼ * Nⱼ₂) - (g * Nⱼ₂)
  du[5] = dNₐ₂ = (g * Nⱼ₂) - (dₐ * Nₐ₂) - (E₂ * Nₐ₂)
  # frac = ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂)
  # du[6] = dS₂ = size_growth_rate * ((S₂ * (1 - ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂))) + (S₁ * ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂))) * (1 - ((S₂ * (1 - ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂))) + (S₁ * ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂))) / (sizeₘₐₓ₂ - (sizeₘₐₓ₂ * E₂)))
  du[6] = dS₂ = size_growth_rate * S₂ * (1 - S₂ / (sizeₘₐₓ₂ - (sizeₘₐₓ₂ * E₂)))
end

p_2 = [0.6, 0.1, 0.05, 0.08, 0.1,
  0.01, 0.04,
  0.05, 0.4,
  1e4, 1e4,
  55.0, 55.0
]
u0_2 = [1e3, 1e3, 40.0,
  1e3, 1e3, 40.0
]
tspan_2 = (0.0, 100.0)
prob_2 = ODEProblem(two_sites!, u0_2, tspan_2, p_2)
sol_2 = solve(prob_2, Tsit5())

plot(sol_2, vars=(0, 1), label="Nⱼ₁")
plot!(sol_2, vars=(0, 2), label="Nₐ₁")
plot!(sol_2, vars=(0, 4), label="Nⱼ₂")
plot!(sol_2, vars=(0, 5), label="Nₐ₂")

plot(sol_2, vars=(0, 3), label="S₁")
plot!(sol_2, vars=(0, 6), label="S₂")


### ----------------------------------------------------------------
### 2.1 General model
### ----------------------------------------------------------------

## Starting the model
#---------------------------
nsites = 10
nlifestages = 5

## initial state of the system
# Initial population sizes at each life stage per site
# There are 100 individuals per m2. We assume 2 km2 per site -> 20k individuals per site.
u0_general = [[18_000.0 for j in 1:nlifestages] for i in 1:nsites]
# Initial average sizes per site. Avg. size is between 45 to 51.
for i in 1:nsites
  push!(u0_general[i], 48.0)
end

# u0 cannot be a nested vector.
u0_general = reduce(vcat, u0_general)

## defining model parameters
function exploit(t, rate)
  if (t % 365) / 365 < 0.42
    return rate
  else
    return 0.0
  end
end

"""
captures the reproductive cycle and equals to zero when fishing and to one when organisms reproduce
"""
function reproductive_cycle(t)
  if (t % 365) / 365 > 0.42
    return 1.0
  else
    return 0.0
  end
end

# conversion rates between stages
# average oocytes per year per adult
avg_oocytes = mean([92098, 804183])
reggs = avg_oocytes / 365  # conversion rate of adults to eggs
r = [reggs, 0.998611, 0.971057, 0.683772, 0.00629]
# natural death rates per life stage.
d = [0.005, 0.005, 0.005, 0.005, 0.000322]
size_growth_rate = 0.32/365
distance_matrix = rand(0.01:0.01:0.1, nsites, nsites)  # TODO use empirical values
distance_matrix[diagind(distance_matrix)] .= 0.0
exploitation_rates = rand(0.001:0.001:0.003, nsites)  # TODO: use empirical values
size_max = 56.0
K = 24_000  # for 2.4 km2 per site. NB: we use K only for the Juveniles.
p_general = [r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K]

function nsites!(du, u, p, t)
  # Change these parameters if you change the model
  nsites = 10
  nlifestages = 5
  site_indices = [0, 6, 12, 18, 24, 30, 36, 42, 48, 54]

  counter = 1
  for site in 1:nsites
    # site_index = (site * (nlifestages + 1)) - (nlifestages + 1)
    site_index = site_indices[site]

    # Stages Egg, Trochophore (1, 2)
    for stage in 1:2
      prev_stage = stage - 1 > 0 ? stage - 1 : nlifestages
      du[counter] = (reproductive_cycle(t) * p[1][stage] * u[site_index+prev_stage]) -
                    (p[1][stage+1] * u[site_index+stage]) -
                    (p[2][stage] * u[site_index+stage]) +
                    (sum(p[4][:, site] .* u[site_indices.+stage])) -
                    (u[site_index+stage]*sum(p[4][site, :]))
      counter += 1
    end

    #stage 3 Veliger
    stage = 3
    prev_stage = 2
    du[counter] = (reproductive_cycle(t) * p[1][stage] * u[site_index+prev_stage]) -
              (p[1][stage+1] * u[site_index+stage]) -
              (p[2][stage] * u[site_index+stage]) +
              (sum(p[4][:, site] .* u[site_indices.+stage])) -
              (u[site_index+stage]*sum(p[4][site, :]))
      counter += 1
      
    # stage 4 Juvenile
    stage = 4
    prev_stage = 3
    du[counter] = (p[1][stage] * u[site_index+prev_stage] * ((p[7] - u[site_index+prev_stage]) / p[7])) -
                (p[1][stage+1] * u[site_index+stage]) - 
                (p[2][stage] * u[site_index+stage]) 
    counter += 1

    # stage 5 adult
    stage = 5
    prev_stage = 4
    du[counter] = (p[1][stage] * u[site_index+prev_stage]) - 
                (exploit(t, p[5][site]) * u[site_index+stage]) -
                (p[2][stage] * u[site_index+stage])
    counter += 1

    # adult sizes
    # adult sizes. dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
    du[counter] = p[3] * u[site_index+nlifestages+1] * (1 - u[site_index+nlifestages+1] / (p[6] - (p[6] * exploit(t, p[4][site]))))
    counter += 1

  end
end

tspan_general = (0.0, 250)
prob_general = ODEProblem(nsites!, u0_general, tspan_general, p_general)
sol_general = solve(prob_general, Tsit5())

plot(sol_general, legend=false)

### --------------------------------------------------
### 3. Sensitivity analysis
### --------------------------------------------------

### 3.1 Single site

f1 = function (p)
  prob1 = remake(prob_1; p=p)
  sol = solve(prob1, Tsit5(); saveat=tspan_1)
  [sol[1, end], sol[2, end], sol[3, end]]
end

bounds_1 = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [1e4, 1e5], [0.0, 1.0], [55.0, 55.0]]

m = gsa(f1, Sobol(), bounds_1, N=10000)

param_names = ["r", "g", "dⱼ", "dₐ", "E", "K", "size_growth_rate", "sizeₘₐₓ"]

p1 = bar(param_names, m.S1[1, :], title="First order Nj", legend=false, xrotation=60)
p2 = bar(param_names, m.ST[1, :], title="Total order Nj", legend=false, xrotation=60)
p3 = bar(param_names, m.S1[2, :], title="First order Na", legend=false, xrotation=60)
p4 = bar(param_names, m.ST[2, :], title="Total order Na", legend=false, xrotation=60)
p5 = bar(param_names, m.S1[3, :], title="First order S", legend=false, xrotation=60)
p6 = bar(param_names, m.ST[3, :], title="Total order S", legend=false, xrotation=60)

plot(p1, p3, p2, p4) # juvenile and adult
plot(p5, p6) # body size

### 3.2 Two sites 

f2 = function (p)
  prob1 = remake(prob_2; p=p)
  sol = solve(prob1, Tsit5(); saveat=tspan_2)
  [sol[1, end], sol[2, end], sol[4, end], sol[5, end]]
end

bounds_2 = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [1e4, 1e5], [1e4, 1e5], [55.0, 55.0], [55.0, 55.0]]

m = gsa(f2, Sobol(), bounds_2, N=10000)

param_names = ["r", "g", "dⱼ", "dₐ", "size_growth_rate", "mₒᵤₜ₁", "mₒᵤₜ₂", "E₁", "E₂", "K₁", "K₂", "sizeₘₐₓ₁", "sizeₘₐₓ₂"]

p1 = bar(param_names, m.S1[1, :], title="First order Nj1", legend=false, xrotation=60)
p2 = bar(param_names, m.ST[1, :], title="Total order Nj1", legend=false, xrotation=60)
p3 = bar(param_names, m.S1[2, :], title="First order Na1", legend=false, xrotation=60)
p4 = bar(param_names, m.ST[2, :], title="Total order Na1", legend=false, xrotation=60)
p5 = bar(param_names, m.S1[3, :], title="First order Nj2", legend=false, xrotation=60)
p6 = bar(param_names, m.ST[3, :], title="Total order Nj2", legend=false, xrotation=60)
p7 = bar(param_names, m.S1[4, :], title="First order Na2", legend=false, xrotation=60)
p8 = bar(param_names, m.ST[4, :], title="Total order Na2", legend=false, xrotation=60)

plot(p1, p3, p2, p4) # Pop 1, juvenile and adult
plot(p5, p7, p6, p8) # Pop 2, juvenile and adult
plot(p1, p2, p5, p6) # Pop 1 vs pop 2, juvenile
plot(p3, p4, p7, p8) # Pop 1 vs pop 2, adults

### -------------------------------
### 4. Data fitting
### --------------------------------

df = CSV.read("data.csv", DataFrame)

sitedfs = groupby(df, "sampling_site")
# sitedfs.keymap

size_by_year_all_sites = DataFrame[]
for site in 1:length(sitedfs)
  site1year = groupby(sitedfs[site], [:year, :species])
  # select(site1year, :total_length_mm)
  cc = combine(site1year, :total_length_mm .=> [mean, std, median])
  push!(size_by_year_all_sites, cc)
end

size.(size_by_year_all_sites, 1)
# Many sites are sampled at one or a few years only.
# Estimate exploitation E for each site only according to one time point.

dfy = groupby(df, [:sampling_site])
species_size_range = combine(dfy, :total_length_mm .=> [maximum, minimum, mean, median, std])
sizeₘₐₓ = maximum(species_size_range.total_length_mm_maximum)
sizeₘᵢₙ = minimum(species_size_range.total_length_mm_minimum)
initial_size = species_size_range.total_length_mm_median
size_growth_rate = 0.2  # It has minimal effect

function estimate_E_across_sites(initial_sizes, sizeₘₐₓ, size_growth_rate)
  sizeonly!(size, E, t) = size_growth_rate * size * (1 - size / (sizeₘₐₓ - (sizeₘₐₓ * E)))
  estimated_Es = zeros(length(initial_sizes))
  for site in 1:length(initial_sizes)
    u0 = initial_sizes[site]
    tspan = (0.0, 60.0)
    p = 0.4872
    prob = ODEProblem(sizeonly!, u0, tspan, p)

    # sol = solve(prob,Tsit5())

    # data = fill(initial_size[1], 60)
    # t = 1:60
    # cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data), maxiters=10000,verbose=false)
    # result = optimize(cost_function, 0.0, 1.0)

    loss_func(sol) = (initial_sizes[site] - sol[end])^2
    cost_function = build_loss_objective(prob, Tsit5(), loss_func, maxiters=10000, verbose=false)
    result = optimize(cost_function, 0.0, 1.0)
    estimated_Es[site] = result.minimizer
  end
  return estimated_Es
end

estimated_E = estimate_E_across_sites(initial_size, sizeₘₐₓ, size_growth_rate)

histogram(estimated_E, bins=30, legend=false, color="black", xlabel="Exploitation", ylabel="Number of Sites")
Plots.savefig("dist_E.pdf")

### -------------------------------
### 0. Testing functions
### --------------------------------
x = 0:0.005:1
y1(x) = x * (1 - x)
y2(x) = x^(2 / 3) - x
y3(x) = -x * log(x)
plot(x, y1.(x), label="Logistic", xlabel="Size", ylabel="Growth rate")
plot!(x, y2.(x), label="von Bertalanffy")
plot!(x, y3.(x), label="Gompertz")
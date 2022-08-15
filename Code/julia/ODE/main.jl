using Pkg
Pkg.activate(".")
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
# In the ODEs, the unit of right hand side is individuals per day. This means that each time step is one day. The rate of the left-hand side are calculated such that after after the given time between life stages, 99% of the individuals from one stage are converted to the next stage. To that end, I solve the followin equation: 0.01 = r^t, where I replace t with the time between stages.
## TODO: Create another version of model where there are five life stages instead of two. In this system, the probability to migrate decreases exponentially. $1/d \times e^{-\gamma}$.
## [x]: Check how to add all equations in a loop instead of writing by hand.

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
### 2.1 N-sites in a loop TODO
### ----------------------------------------------------------------


# dummy ode to test (It works!)
function nsites!(du, u, p, t)
  for site in 1:5
    du[site] = p[site] * u[site]
  end
end

pn = [0.1, 0.2, 0.3, 0.4, 0.5]
u0n = [0.15, 0.14, 0.13, 0.12, 0.11]
tspann = (0.0, 10.0)
probn = ODEProblem(nsites!, u0n, tspann, pn)
soln = solve(probn, Tsit5())

plot(soln)

# Real ODE
function nsites!(du, u, p, t)

end


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
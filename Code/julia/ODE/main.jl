using Pkg
Pkg.activate(".")
using OrdinaryDiffEq
# using DifferentialEquations
using GlobalSensitivity
using Plots
using Statistics
using DataFrames
using CSV

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
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E)))
end

p_1 = [0.6, 0.6, 0.05, 0.08, 0.2, 1e4, 0.2, 40.0]
u0_1 = [1e3, 1e3, 40.0]
tspan_1 = (0.0, 50.0)
prob_1 = ODEProblem(single_site!, u0_1, tspan_1, p_1)
sol_1 = solve(prob_1, Tsit5())

plot(sol_1, vars = (0, 1), label = "Nⱼ")
plot!(sol_1, vars = (0, 2), label = "Nₐ")

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

plot(sol_2, vars = (0, 1), label = "Nⱼ₁")
plot!(sol_2, vars = (0, 2), label = "Nₐ₁")
plot!(sol_2, vars = (0, 4), label = "Nⱼ₂")
plot!(sol_2, vars = (0, 5), label = "Nₐ₂")

plot(sol_2, vars = (0, 3), label = "S₁")
plot!(sol_2, vars = (0, 6), label = "S₂")

### --------------------------------------------------
### 3. Sensitivity analysis
### --------------------------------------------------

### 3.1 Single site

f1 = function (p)
  prob1 = remake(prob_1; p = p)
  sol = solve(prob1, Tsit5(); saveat = tspan_1)
  [sol[1, end], sol[2, end], sol[3, end]]
end

bounds_1 = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [1e4, 1e5], [0.0, 1.0], [55.0, 55.0]]

m = gsa(f1, Sobol(), bounds_1, N = 10000)

param_names = ["r", "g", "dⱼ", "dₐ", "E", "K", "size_growth_rate", "sizeₘₐₓ"]

p1 = bar(param_names, m.S1[1, :], title = "First order Nj", legend = false, xrotation = 60)
p2 = bar(param_names, m.ST[1, :], title = "Total order Nj", legend = false, xrotation = 60)
p3 = bar(param_names, m.S1[2, :], title = "First order Na", legend = false, xrotation = 60)
p4 = bar(param_names, m.ST[2, :], title = "Total order Na", legend = false, xrotation = 60)
p5 = bar(param_names, m.S1[3, :], title = "First order S", legend = false, xrotation = 60)
p6 = bar(param_names, m.ST[3, :], title = "Total order S", legend = false, xrotation = 60)

plot(p1, p3, p2, p4) # juvenile and adult
plot(p5,p6) # body size

### 3.2 Two sites 

f2 = function (p)
  prob1 = remake(prob_2; p = p)
  sol = solve(prob1, Tsit5(); saveat = tspan_2)
  [sol[1, end], sol[2, end], sol[4, end], sol[5, end]]
end

bounds_2 = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [1e4, 1e5], [1e4, 1e5], [55.0, 55.0], [55.0, 55.0]]

m = gsa(f2, Sobol(), bounds_2, N = 10000)

param_names = ["r", "g", "dⱼ", "dₐ", "size_growth_rate", "mₒᵤₜ₁", "mₒᵤₜ₂", "E₁", "E₂", "K₁", "K₂", "sizeₘₐₓ₁", "sizeₘₐₓ₂"]

p1 = bar(param_names, m.S1[1, :], title = "First order Nj1", legend = false, xrotation = 60)
p2 = bar(param_names, m.ST[1, :], title = "Total order Nj1", legend = false, xrotation = 60)
p3 = bar(param_names, m.S1[2, :], title = "First order Na1", legend = false, xrotation = 60)
p4 = bar(param_names, m.ST[2, :], title = "Total order Na1", legend = false, xrotation = 60)
p5 = bar(param_names, m.S1[3, :], title = "First order Nj2", legend = false, xrotation = 60)
p6 = bar(param_names, m.ST[3, :], title = "Total order Nj2", legend = false, xrotation = 60)
p7 = bar(param_names, m.S1[4, :], title = "First order Na2", legend = false, xrotation = 60)
p8 = bar(param_names, m.ST[4, :], title = "Total order Na2", legend = false, xrotation = 60)

plot(p1, p3, p2, p4) # Pop 1, juvenile and adult
plot(p5, p7, p6, p8) # Pop 2, juvenile and adult
plot(p1, p2, p5, p6) # Pop 1 vs pop 2, juvenile
plot(p3, p4, p7, p8) # Pop 1 vs pop 2, adults

### -------------------------------
### 4. Data fitting
### --------------------------------

df = CSV.read("data.csv", DataFrame)

sitedfs = groupby(df, "sampling_site")
sitedfs.keymap

size_by_year_all_sites = DataFrame[]
for site in 1:length(sitedfs)
  site1year = groupby(sitedfs[site], [:year, :species]);
  # select(site1year, :total_length_mm)
  cc = combine(site1year, :total_length_mm .=> [mean,std, median])
  push!(size_by_year_all_sites, cc)
end

size.(size_by_year_all_sites, 1)
# Many sites are sampled at one or a few years only.

function sizeonly!(du, u, p, t)
  Sₐ = u
  size_growth_rate, sizeₘₐₓ, E = p
  du[1] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E)))
end


### -------------------------------
### 0. Testing functions
### --------------------------------
x = 0:0.005:1
y1(x) = x * (1 - x)
y2(x) = x^(2 / 3) - x
y3(x) = -x * log(x)
plot(x, y1.(x), label = "Logistic", xlabel = "Size", ylabel = "Growth rate")
plot!(x, y2.(x), label = "von Bertalanffy")
plot!(x, y3.(x), label = "Gompertz")
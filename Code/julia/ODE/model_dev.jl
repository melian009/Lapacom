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

## [x]: Empirically estimate the rate at which individuals from one stage turn into another one:
## Conversion rates for the given number of days is as follows:
# 730: 0.00629
# 3: 0.784557
# 0.7: 0.998611
# 1.3: 0.971057
# 7: 0.4820525
# In the ODEs, the unit of left side is individuals per day. This means that each time step is one day. The rate of the right side are calculated such that after the given time between life stages, 99% of the individuals from one stage are converted to the next stage. To that end, I solve the following equation: 0.01 = r^t, where I replace t with the time between stages.
# [x]: Modify the equations according to the following facts:
# 1. adults live 7 to 10 years, depending on their size. 
#   Age in days ranges between 2555 to 3650.
#   Actual age is (size * max age) / max size. If actual age < min age, actual age = min age. This is difficult to include in the ODE. We can ignore it and use an average age: 3102 days.
#   death probability of adults da =  1 / age. so 1/3102 = 0.000322.
# 2. adults spawn 92,098 to 804,183 oocytes per year during the spawning season.
#   Appearance rate of eggs: #oocytes/365 in a sin function (to account for spawning season). Note that adults do not turn into eggs (they live many steps.)
## [x]: estimate `size_growth_rate` by knowing juvenile average size and adult average size and dividing their difference by 730 days. Limpets 0.32 cm per year, which is 0.32/365 cm per day. To calculate the size_growth_rate, see estimate_size_growth_rate.jl.
## [x] Estimate exploitation rates given average sizes and a single pool model.



#=
Five simplifying assumptions:

1. Transition times. Some short transition times can be aggregated.
2. The number of oocytes. We use the mean number of oocytes laid per adult per year.
3. Dispersal potential of trochophore and veliger. Eggs disperse passively, T and V can disperse actively.
4. Juveniles homing. Eggs have a global pool. Once they are released, they travel with the currents to the ocean. In the ocean they develop and actively migrate back to different sites. This is like the "lottery model".
5. Exploitation. Seven months of exploitation (April to October), five months rest (November to March). There is no reproduction during the exploitation time. NOTE that adults continueously keep growing.
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

Et(t) = (sin(t * 0.1)^2) / 2  # time varying exploitation
# To adjust the sin function horizonally (stretch/shrink the wave length), multiply t by a factor c.

# # You may also use a function with if statements.
# function Et(t)
#   if modf(t)[1] < 0.5
#     return 0.0
#   else
#     return 0.5
#   end
# end

# # Using this function results in errors.
# function Et(t)
#   if (t % 365) / 365 >= 0.42
#     return 1.0
#   else
#     return 5.0
#   end
# end

p_1 = [5.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0]
u0_1 = [1e3, 1e3, 40.0]
tspan_1 = (0.0, 500.0)
prob_1 = ODEProblem(single_site!, u0_1, tspan_1, p_1)
sol_1 = solve(prob_1, Tsit5())

fig, ax, plt = lines(sol_1.t, [i[1] for i in sol_1.u], label="Nⱼ")
lines!(ax, sol_1.t, [i[2] for i in sol_1.u], label="Nₐ")
fig[1, 2] = Legend(fig, ax, "Stage")
save("one_site_model.png", fig)

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
  du[2] = dNₐ₁ = (g * Nⱼ₁) - (dₐ * Nₐ₁) - (E₁(t) * Nₐ₁)
  # track the influence of migration on body size. If there is immigration from a site where individuals are genetically small, then it should affect the distribution of body size in this site too. 
  # S₁ need to be replaced by the mean body size of this population which is a function of body sizes of all source populations and their contribution to this sink population.
  # frac = ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁)
  # du[3] = dS₁ = size_growth_rate * ((S₁ * (1 - ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁))) + (S₂ * ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁))) * (1 - ((S₁ * (1 - ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁))) + (S₂ * ((mₒᵤₜ₂ * Nⱼ₂) / Nⱼ₁))) / (sizeₘₐₓ₁ - (sizeₘₐₓ₁ * E₁)))
  du[3] = dS₁ = size_growth_rate * S₁ * (1 - S₁ / (sizeₘₐₓ₁ - (sizeₘₐₓ₁ * E₁(t))))  # logistic growth: rN(1- N/K)

  du[4] = dNⱼ₂ = (r * Nₐ₂ * ((K₂ - Nₐ₂) / K₂)) + (mₒᵤₜ₁ * Nⱼ₁) - (mₒᵤₜ₂ * Nⱼ₂) - (dⱼ * Nⱼ₂) - (g * Nⱼ₂)
  du[5] = dNₐ₂ = (g * Nⱼ₂) - (dₐ * Nₐ₂) - (E₂(t) * Nₐ₂)
  # frac = ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂)
  # du[6] = dS₂ = size_growth_rate * ((S₂ * (1 - ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂))) + (S₁ * ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂))) * (1 - ((S₂ * (1 - ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂))) + (S₁ * ((mₒᵤₜ₁ * Nⱼ₁) / Nⱼ₂))) / (sizeₘₐₓ₂ - (sizeₘₐₓ₂ * E₂)))
  du[6] = dS₂ = size_growth_rate * S₂ * (1 - S₂ / (sizeₘₐₓ₂ - (sizeₘₐₓ₂ * E₂(t))))
end

Et1(t) = (sin(t * 0.1)^2) / 1.1
Et2(t) = (sin(t * 0.1)^2) / 2

p_2 = [6.6, 0.1, 0.05, 0.08, 0.1,
  0.01, 0.04,
  Et1, Et2,
  1e4, 1e4,
  55.0, 55.0
]
u0_2 = [1e3, 1e3, 40.0,
  1e3, 1e3, 40.0
]
tspan_2 = (0.0, 100.0)
prob_2 = ODEProblem(two_sites!, u0_2, tspan_2, p_2)
sol_2 = solve(prob_2, Tsit5())


fig, ax, plt = lines(sol_2.t, [i[1] for i in sol_2.u], label="Nⱼ₁")
lines!(ax, sol_2.t, [i[2] for i in sol_2.u], label="Nₐ₁")
lines!(ax, sol_2.t, [i[4] for i in sol_2.u], label="Nⱼ₂")
lines!(ax, sol_2.t, [i[5] for i in sol_2.u], label="Nₐ₂")
fig[1, 2] = Legend(fig, ax, "Stage")
save("two_sites_model.png", fig)


fig, ax, plt = lines(sol_2.t, [i[3] for i in sol_2.u], label="S₁")
lines!(ax, sol_2.t, [i[6] for i in sol_2.u], label="S₂")
fig[1, 2] = Legend(fig, ax, "Species size")
save("two_sites_model_sizes.png", fig)



### ----------------------------------------------------------------
### 2.1 General model
### ----------------------------------------------------------------

datadir = "./"

## Starting the model
#---------------------------
nsites = 8
nlifestages = 5

## initial state of the system
# Initial population sizes at each life stage per site
# There are 100 individuals per m2. We assume 2 km2 per site -> 20k individuals per site.
u0_general = [[58_000.0 for j in 1:nlifestages] for i in 1:nsites]
# Initial average sizes per site. Avg. size is between 45 to 51.
for i in 1:nsites
  push!(u0_general[i], 48.0)  # TODO: how to initialize the avg sizes?
end
# u0 cannot be a nested vector.
u0_general = reshape(reduce(vcat, u0_general), length(u0_general[1]), length(u0_general))'

## defining model parameters

"""
captures the reproductive cycle and equals to zero when fishing and to one when organisms reproduce. Reproduction occurs between November and March (0.42 of the year).
"""
function reproductive_cycle(t)
  if (t % 365) / 365 >= 0.42
    return 0.0
  else
    return 1.0
  end
end

# conversion rates between stages
# average oocytes per year per adult
avg_oocytes = 385_613 # This is the actual mean. mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r = [reggs, 0.998611, 0.971057, 0.4820525, 0.00629]
# natural death rates per life stage.
d = [0.999 / 365, 0.585 / 365, 0.343 / 365, 0.201 / 365, 0.000322]  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = 0.00014898749263737575

distance_df = CSV.read(joinpath(datadir, "distance_matrix.csv"), DataFrame)
distance_matrix = Float64.(Matrix(distance_df)[:, 1:end-1])
distance_matrix = distance_matrix[1:nsites, 1:nsites]
migraation_probs_df = CSV.read(joinpath(datadir, "migration_probabilites_among_sites.csv"), DataFrame)
# Modify distances (migration probs) by movement prob. which are rough estimates based on the land shape and oceanic currents
mig_probs = Float64.(Matrix(migraation_probs_df)[:, 1:end-1])
mig_probs = mig_probs[1:nsites, 1:nsites]


# exploitation_rates = [0.436795998061044, 0.43767209156532155, 0.4603254055329175, 0.38748435327632225, 0.40337922500828566, 0.5105131417482706, 0.4799123913754184, 0.47959950031955256]  # These values are estimated based on the average size of the limpets at each site. See the text below.
exploitation_rates = [0.9880685477383748, 0.9880947574980841, 0.9891331673866427, 0.9874718910553065, 0.9875900832783625, 0.9900632400512024, 0.9894700088322096, 0.9894645471975891] # These values are estimated based on the average size of the limpets at each site. See the text below. These second exploitation rates used Gompertz growth function instead of logistic.
exploitation_rates = exploitation_rates .- exploitation_rates[4]  # subtract the exploitation rate of Desertas from all sites because it is a fully protected area and no fishing happens there. This is to make sure that the exploitation rate of Desertas is zero and the exploitation rates of other sites are relative to Desertas.

# function exploitation_rate(size, site)
#   size_max = 56.0
#   min_size_for_fishing = 40
#   exploitation_rates_max = [0.43679599501684097, 0.43767209022758713, 0.46032540377965536, 0.38748435785424035, 0.4033792240083418, 0.5105131408815863, 0.47991239051754064, 0.4795994994034682]  # These values are estimated based on the average size of the limpets at each site.
#   # exploitation_rates_max = [0.980700696804078, 0.9807285958927162, 0.989105147959682, 0.9796381668584575, 0.9875759368524103, 0.9900648704172368, 0.9894851154977173, 0.9894795735668236] # These exploitation rates used Gompertz growth function instead of logistic.
#   exploitation_rates_max = exploitation_rates_max .- exploitation_rates_max[4]  # subtract the exploitation rate of Desertas from all sites because it is a fully protected area and no fishing happens there. This is to make sure that the exploitation rate of Desertas is zero and the exploitation rates of other sites are relative to Desertas.
#   if size <= min_size_for_fishing
#     return 0.0
#   else
#     slope = 1 / (size_max - min_size_for_fishing)
#     intercept = -min_size_for_fishing / (size_max - min_size_for_fishing)
#     E = slope * size + intercept
#     E = min(E, 1.0)
#     return E * exploitation_rates_max[site]
#   end
# end

size_max = 56.0
K = 64_000  # for 6.4 km2 per site.
α = [0.1, 0.1, 0.1]  # dispersion factor for Egg, Trochophore, and Veliger
# Since we do not have any info about site size, dispersion is only a function of dispersion factor and distance.
p_general = [r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α, mig_probs]

"""
I flatten the parameters and then return them back in the original shape inside `nsites!`. This is because to run Sensitivity analysis (`gsa`), I have to pass flat bounds list.
"""
function flatten_p(p)
  flat_p = Real[]
  for p_list in p
    if typeof(p_list) <: Real
      push!(flat_p, p_list)
    else
      for p in p_list
        push!(flat_p, p)
      end
    end
  end
  return flat_p
end

p_general_flat = flatten_p(p_general)

function restructure_flat_p(flat_params, original_shape)
  nested_params = []
  start_idx = 1
  for p_list in original_shape
    if typeof(p_list) <: Real
      push!(nested_params, flat_params[start_idx])
      start_idx += 1
    else
      end_idx = start_idx + length(p_list) - 1
      push!(nested_params, reshape(eltype(p_list).(flat_params[start_idx:end_idx]), size(p_list)))
      start_idx = end_idx + 1
    end
  end
  return nested_params
end

"""
  Estimates the `size at first maturity` by constructing a linear line between `maximum size` and average size before and after protection (These values come from empirical data).

  Reproduction capacity depends on size where at maximum size, reproduction is 100% and at size at first maturity it is 50%.

  NOTE The function only implements the calculations for P. ordinaria. Write a new function fo P. aspera.

  ## Data


  ### Size at first maturity

  Size at first Maturity bef/after [P. ordinaria,P. aspera]  = [33.4/37.4, 34.6/37.5]

  ### Average size

  Before (1996-2006): FULL ACCESS.  # Use this for before
  Patella apera = 43.53mm
  Patella ordinaria = 46.26mm

  Only MPA  # Use this for after
  After (2007-2017)
  Patella aspera = 50.61mm
  Patella ordinaria = 49.25mm

  ## The model

  Let's represent size at first maturity by "M" and the average body size by "A". To establish the relationship between the two sizes, we can use a simple linear equation:

  M = kA + b

  where "k" is the proportionality constant and "b" is a constant term.

  We have two sets of values, one for the time before fishing protection (M₁/A₁ = 33.4/46.26) and one for the time after fishing protection (M₂/A₂ = 37.4/49.25). To find the values of "k" and "b", we can set up two equations:

  33.4 = k * 46.26 + b
  37.4 = k * 49.25 + b
  Solving this system of linear equations gives us values of "k" and "b".

  M ≈ 1.34A - 28.06
"""
function calculate_size_at_first_maturity(current_avg_size)
  M = 1.34 * (current_avg_size) - 28.06
end

"Return the reproduction capacity (between 0 and 1) given the current average size and size at first maturity and maximum size"
reproduction_capacity(Saverage, Smaturity, Smax) = min(max(0.5 * (1.0 + (Saverage - Smaturity) / (Smax - Smaturity)), 0.0), 1.0)

function nsites!(du, u, p, t)
  original_shape = Any[[5.846581865622962, 0.998611, 0.971057, 0.4820525, 0.00629], [0.001, 0.001, 0.001, 0.001, 0.000322], 0.0008767123287671233, [0.0 12.842542485919632 33.20921059959618 73.75112848881275 45.78816486689583 40.95386826936409 23.127098263365838 13.53847480618783; 12.842542485919632 0.0 31.551937191782326 74.4358482251258 50.57366923475108 43.525278937384996 18.68655292046152 19.67083747567315; 33.20921059959618 31.551937191782326 0.0 43.012716906931786 25.586975451710284 15.767470469506062 13.16139561432816 21.70783981938122; 73.75112848881275 74.43584822512581 43.012716906931786 0.0 32.56225153468865 33.114636105203594 56.164462942582226 60.38662955616323; 45.78816486689583 50.573669234751094 25.586975451710284 32.56225153468865 0.0 10.08940978792187 35.940503401029886 32.44536267645653; 40.95386826936408 43.525278937384996 15.767470469506062 33.114636105203594 10.08940978792187 0.0 27.21061174752855 27.454025172124666; 23.127098263365838 18.68655292046152 13.16139561432816 56.164462942582226 35.940503401029886 27.21061174752855 0.0 16.029653650497625; 13.53847480618783 19.67083747567315 21.70783981938122 60.38662955616323 32.44536267645653 27.454025172124666 16.029653650497625 0.0], [0.002, 0.001, 0.002, 0.001, 0.001, 0.002, 0.002, 0.002], 56.0, 64000, [0.1, 0.1, 0.1], [0.0 0.8 0.0 0.0 0.0 0.0 0.2 0.8; 0.0 0.0 0.2 0.0 0.0 0.2 0.5 0.0; 0.0 0.6 0.0 0.4 0.2 0.4 0.8 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.2 0.8 0.8 0.0 1.0 0.6 0.0; 0.0 0.4 0.8 0.6 0.2 0.0 0.6 0.0; 0.0 0.8 0.4 0.2 0.2 0.4 0.0 0.0; 0.8 0.6 0.0 0.0 0.0 0.0 0.0 0.0]]
  p = restructure_flat_p(p, original_shape)
  nsites = 8
  nlifestages = 5
  reduce_trochophore_dispersals = 0.2 # dispersal rates of Trochophores should be a fraction of Eggs. This factor is by how much it is reduced.

  r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α, mig_probs = p
  for site in 1:nsites

    # Stages 1 Egg
    stage = 1
    prev_stage = 5

    # immigration
    dispersal_probs1 = (exp(-α[stage]) ./ distance_matrix[:, site])
    dispersal_probs1[site] = 0.0
    dispersal_probs1 = dispersal_probs1 .* mig_probs[:, site]
    # dispersal_probs1 = dispersal_probs1 ./ sum(dispersal_probs1)

    # emigration
    dispersal_probs2 = (exp(-α[stage]) ./ distance_matrix[site, :])
    dispersal_probs2[site] = 0.0
    dispersal_probs2 = dispersal_probs2 .* mig_probs[site, :]
    # dispersal_probs2 = dispersal_probs2 ./ sum(dispersal_probs2)

    Saverage = du[site, 6]  # 6 is the index of average size
    Smaturity = calculate_size_at_first_maturity(Saverage)

    du[site, stage] = (reproductive_cycle(t) * r[stage] * u[site, prev_stage] * reproduction_capacity(Saverage, Smaturity, size_max)) -  
                      (r[stage+1] * u[site, stage]) -
                      (d[stage] * u[site, stage]) +
                      (sum(dispersal_probs1 .* u[:, stage])) -
                      (u[site, stage] * sum(dispersal_probs2))
    # What does the last part of the first term do? We use it to reduce the reproductive capacity of the adults based on their size. At max size, 100% of them reproduce. At "size at first maturity" size, only 50% of them reproduce. Make a ramp function that calculates the fraction of the adults that reproduce given the mean size of the individuals. For this, we will also need to have another equation (`calculate_size_at_first_maturity`) that calculates the "size at the first maturity" because its a function of the mean size of the individuals. For that, we use a linear relationship between mean size and size at the first maturity before/after introducing fishing protection.

    # Stage 2 Trochophore
    stage = 2
    prev_stage = 1

    # immigration
    dispersal_probs1 = (exp(-α[stage]) ./ distance_matrix[:, site])
    dispersal_probs1[site] = 0.0
    dispersal_probs1 = dispersal_probs1 .* mig_probs[:, site] .* reduce_trochophore_dispersals
    # dispersal_probs1 = dispersal_probs1 ./ sum(dispersal_probs1)

    # emigration
    dispersal_probs2 = (exp(-α[stage]) ./ distance_matrix[site, :])
    dispersal_probs2[site] = 0.0
    dispersal_probs2 = dispersal_probs2 .* mig_probs[site, :] .* reduce_trochophore_dispersals
    # dispersal_probs2 = dispersal_probs2 ./ sum(dispersal_probs2)

    du[site, stage] = (r[stage] * u[site, prev_stage]) -
                      (r[stage+1] * u[site, stage]) -
                      (d[stage] * u[site, stage]) +
                      (sum(dispersal_probs1 .* u[:, stage])) -
                      (u[site, stage] * sum(dispersal_probs2))

    #stage 3 Veliger
    stage = 3
    prev_stage = 2

    du[site, stage] = (r[stage] * u[site, prev_stage] * ((K - u[site, stage]) / K)) -
                      (r[stage+1] * u[site, stage]) -
                      (d[stage] * u[site, stage])

    # stage 4 Juvenile
    stage = 4
    prev_stage = 3
    du[site, stage] = (r[stage] * u[site, prev_stage] * ((K - u[site, stage]) / K)) -
                      (r[stage+1] * u[site, stage]) -
                      (d[stage] * u[site, stage])

    # stage 5 adult
    stage = 5
    prev_stage = 4
    du[site, stage] = (r[stage] * u[site, prev_stage] * ((K - u[site, stage]) / K)) -
                      (((1 - reproductive_cycle(t)) * exploitation_rates[site]) * u[site, stage]) -
                      (d[stage] * u[site, stage])

    # adult sizes
    stage = 6
    du[site, stage] = size_growth_rate * u[site, nlifestages+1] * (1 - (u[site, nlifestages+1] / (size_max - (size_max * ((1 - reproductive_cycle(t)) * exploitation_rates[site])))))
    # du[site, stage] = size_growth_rate * u[site, nlifestages+1] * exp(-u[site, nlifestages+1] / (size_max - (size_max * ((1 - reproductive_cycle(t)) * exploitation_rates[site]))))
  end
  # for i in 1:nsites
  #   for j in 1:(nlifestages-1)
  #     du[i, j] = max(du[i, j], 0.0)
  #   end
  # end
end

tspan_general = (0.0, 3000.0)
prob_general = ODEProblem(nsites!, u0_general, tspan_general, p_general_flat)
sol_general = solve(prob_general)

# Plot
all_u = sol_general.u
all_times = sol_general.t
# site_names = distance_df.site
site_names = ["Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente"]
for stage in 1:5
  fig = Figure()
  ax1 = Axis(fig[1, 1])
  lines!(ax1, all_times, [all_u[i][1, stage] for i in 1:length(all_times)], yscale=:log10, label=site_names[1])
  ax1.title = "N for stage: $(stage)"
  for site in 2:nsites#2:8
    lines!(ax1, all_times, [all_u[i][site, stage] for i in 1:length(all_times)], label=site_names[site])
  end
  fig[1, 2] = Legend(fig, ax1, "Site")
  save("figs/stage=$stage.pdf", fig)
end

# Changes in body size
stage = 6
fig, ax, plt = lines(all_times, [all_u[i][1, stage] for i in 1:length(all_times)], yscale=:log10, label=site_names[1])
ax.title = "Body size"
for site in 2:nsites
  lines!(ax, all_times, [all_u[i][site, stage] for i in 1:length(all_times)], label=site_names[site])
end
fig[1, 2] = Legend(fig, ax, "Site")
save("figs/body_sizes.pdf", fig)

### --------------------------------------------------
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

m = gsa(f2, Sobol(), bounds_2, N=100, batch=true)

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

### 3.3 n sites

lifestage = 5

# for factor in e_factors
factor = 1.0

exploitation_rates = exploitation_rates_org * factor
p_general = [r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α, mig_probs]
p_general_flat = flatten_p(p_general)

start_year = 2006
end_year = 2018
nyears = end_year - start_year + 1
ndays = nyears * 365.0
tspan_general = (0.0, ndays)

prob_general = ODEProblem(nsites!, u0_general, tspan_general, p_general_flat)
# sol_general = solve(prob_general)

fn = function (p)
  prob1 = remake(prob_general; p=p)
  sol = solve(prob1)
  indices = (lifestage - 1) * 8 .+ (1:8)  # 8 sites
  [mean(sol[indices[1], :]), mean(sol[indices[2], :]), mean(sol[indices[3], :]), mean(sol[indices[4], :]), mean(sol[indices[5], :]), mean(sol[indices[6], :]), mean(sol[indices[7], :]), mean(sol[indices[8], :])] # mean population size for each site across all times.
end

mig_probs_bounds = Array{Float64}[]
for ind in eachindex(mig_probs)
  entry = [mig_probs[ind], mig_probs[ind] + 0.001]
  push!(mig_probs_bounds, entry)
end
distance_matrix_bounds = Array{Float64}[]
for ind in eachindex(distance_matrix)
  entry = [distance_matrix[ind], distance_matrix[ind] + 0.001]
  push!(distance_matrix_bounds, entry)
end
bounds_n = [
  [0.0, 6], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0],
  [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0],
  [0.0, 0.1],
  distance_matrix_bounds...,
  [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0],
  [16.0, 76.0],
  [20000, 90000],
  [0.0, 0.2], [0.0, 0.2], [0.0, 0.2],
  mig_probs_bounds...
]

output = "../figs/sensitivity_results_lifestage=$(lifestage)_e_factor=$(factor).jld2"
isfile(output) ? load(output, "results") : m = gsa(fn, Sobol(), bounds_n, samples=100, batch=true) && save(output, "results", m)

# Plot the results
for site in 1:8
  param_names = ["r", "d", "size_growth_rate", "distance_matrix", "exploitation_rates", "size_max", "K", "α", "mig_probs"]
  important_parameter_indices = vcat(collect(1:5), # r for all life stages
    vcat(6:10), # d for all life stages
    vcat(76:83), # exploitation rates for all sites. (10+64+1):(10+64+1+8)
    vcat(86:88), # α
  )
  important_parameter_names = ["r1", "r2", "r3", "r4", "r5", "d1", "d2", "d3", "d4", "d5", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "α1", "α2", "α3"]  ## Bar plot
  normalized_mST = m.ST ./ sum(m.ST, dims=2)
  sensitivity_indices = normalized_mST[site, important_parameter_indices]

  f = Figure()
  ax = Axis(f[1, 1], ylabel="Sensitivity index - total effect", xlabel="Parameter", xticks=(1:length(important_parameter_names), important_parameter_names))
  bars = barplot!(ax, 1:length(important_parameter_names), sensitivity_indices, color=:black)
  save("../figs/sensitivity_indices_barplot_total_effect_site=$(site)_lifestage=$(lifestage)_e_factor=$(factor).pdf", f)
end

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
size_growth_rate = 0.00014898749263737575  # It has minimal effect

function estimate_E_across_sites(initial_sizes, sizeₘₐₓ, size_growth_rate)
  sizeonly!(size, E, t) = size_growth_rate * size * (1 - size / (sizeₘₐₓ - (sizeₘₐₓ * E)))
  # sizeonly!(size, E, t) = size_growth_rate * size * exp(-size / (sizeₘₐₓ - (sizeₘₐₓ * E)))
  estimated_Es = zeros(length(initial_sizes))
  for site in 1:length(initial_sizes)
    u0 = initial_sizes[site]
    tspan = (0.0, 60.0)
    p = 0.4872
    prob = ODEProblem(sizeonly!, u0, tspan, p)

    loss_func(sol) = (initial_sizes[site] - sol[end])^2
    cost_function = build_loss_objective(prob, Tsit5(), loss_func, maxiters=10000, verbose=false)
    result = optimize(cost_function, 0.0, 1.0)
    estimated_Es[site] = result.minimizer
  end
  return estimated_Es
end

estimated_E = estimate_E_across_sites(initial_size, sizeₘₐₓ, size_growth_rate)
# For the logistic growth model: [0.43679599501684097, 0.43767209022758713, 0.46032540377965536, 0.38748435785424035, 0.4033792240083418, 0.5105131408815863, 0.47991239051754064, 0.4795994994034682]
# For the Gompertz growth model: [0.980700696804078, 0.9807285958927162, 0.989105147959682, 0.9796381668584575, 0.9875759368524103, 0.9900648704172368, 0.9894851154977173, 0.9894795735668236]
fig = Figure()
ax = Axis(fig)
p = hist!(ax, estimated_E, bins=30, color=:black, xlabel="Exploitation", ylabel="Number of Sites")

save("figs/estimated_exploration_rate_per_site.pdf", p)

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
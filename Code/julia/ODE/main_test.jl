using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim




### ----------------------------------------------------------------
### 2.1 General model
### ----------------------------------------------------------------

## Starting the model
#---------------------------
nsites = 8
nlifestages = 5

## initial state of the system
# Initial population sizes at each life stage per site
# There are 100 individuals per m2. We assume 2 km2 per site -> 20k individuals per site.
u0_general = [[1_800.0 for j in 1:nlifestages] for i in 1:nsites]
# Initial average sizes per site. Avg. size is between 45 to 51.
for i in 1:nsites
  push!(u0_general[i], 48.0)
end
# u0 cannot be a nested vector.
u0_general = reshape(reduce(vcat, u0_general), length(u0_general[1]), length(u0_general))'

## defining model parameters

"""
captures the reproductive cycle and equals to zero when fishing and to one when organisms reproduce
"""
function reproductive_cycle(t)
  if (t % 365) / 365 >= 0.42
    return 1.0
  else
    return 0.0
  end
end

# conversion rates between stages
# average oocytes per year per adult
```Oocyte values ​are for Patella ordinaria only. 
It needs to include the values ​​for the other species.```
avg_oocytes = mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
reggs = reggs / 500 # because the rate is too high to be handled

r = [reggs, 0.998611, 0.971057, 0.4820525, 0.00629]

# natural death rates per life stage.
d = [0.001, 0.001, 0.001, 0.001, 0.000322]
# d = 

#=
Empirical estimated mortality (Z,d,F) and exploitation rates (E):

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

For numerical simulation using these rates d and E.

Mean size: 
Before (1996-2006)
P apera = 43.53mm
P ordinaria = 46.26mm
After  (2007-2017)
P aspera = 44.45mm
P ordinaria = 46.44mm

=#

size_growth_rate = 0.32 / 365
distance_df = CSV.read("distance_matrix.csv", DataFrame)
distance_matrix = Float64.(Matrix(distance_df)[:, 1:end-1])
distance_matrix = distance_matrix[1:nsites, 1:nsites]
exploitation_rates = rand(0.001:0.001:0.003, nsites)  # TODO: use empirical values
size_max = 56.0
K = 64_000  # for 6.4 km2 per site.
α = [0.1, 0.1, 0.1]  # dispersion factor for Egg, Trochophore, and Veliger
# Since we do not have any info about site size, dispersion is only a function of dispersion factor and distance.
p_general = [r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α]

function nsites!(du, u, p, t)
  nsites = 1
  nlifestages = 5

  r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α = p
  for site in 1:nsites

    # Stages 1 Egg
    stage = 1
    prev_stage = 5

    dispersal_probs1 = (exp(-α[stage]) ./ distance_matrix[:, site])
    dispersal_probs1[site] = 0.0
    dispersal_probs2 = (exp(-α[stage]) ./ distance_matrix[site, :])
    dispersal_probs2[site] = 0.0

    du[site, stage] = (reproductive_cycle(t) * r[stage] * u[site, prev_stage]) -
                      (r[stage+1] * u[site, stage]) -
                      (d[stage] * u[site, stage]) +
                      (sum(dispersal_probs1 .* u[:, stage])) -
                      (u[site, stage] * sum(dispersal_probs2))

    # Stage 2 Trochophore
    stage = 2
    prev_stage = 1

    dispersal_probs1 = (exp(-α[stage]) ./ distance_matrix[:, site])
    dispersal_probs1[site] = 0.0
    dispersal_probs2 = (exp(-α[stage]) ./ distance_matrix[site, :])
    dispersal_probs2[site] = 0.0

    du[site, stage] = (r[stage] * u[site, prev_stage]) -
                      (r[stage+1] * u[site, stage]) -
                      (d[stage] * u[site, stage]) +
                      (sum(dispersal_probs1 .* u[:, stage])) -
                      (u[site, stage] * sum(dispersal_probs2))

    #stage 3 Veliger
    stage = 3
    prev_stage = 2

    dispersal_probs1 = (exp(-α[stage]) ./ distance_matrix[:, site])
    dispersal_probs1[site] = 0.0
    dispersal_probs2 = (exp(-α[stage]) ./ distance_matrix[site, :])
    dispersal_probs2[site] = 0.0

    du[site, stage] = (r[stage] * u[site, prev_stage] * ((K - u[site, stage]) / K)) -
                      (r[stage+1] * u[site, stage]) -
                      (d[stage] * u[site, stage]) +
                      (sum(dispersal_probs1 .* u[:, stage])) -
                      (u[site, stage] * sum(dispersal_probs2))

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
                      ((reproductive_cycle(t) * exploitation_rates[site]) * u[site, stage]) -
                      (d[stage] * u[site, stage])

    # adult sizes
    stage = 6
    du[site, stage] = size_growth_rate * u[site, nlifestages+1] * (1 - u[site, nlifestages+1] / (size_max - (size_max * (reproductive_cycle(t) * exploitation_rates[site]))))
  end
end

tspan_general = (1.0, 3000.0)
prob_general = ODEProblem(nsites!, u0_general, tspan_general, p_general)
sol_general = solve(prob_general)

# Plot
all_u = sol_general.u
all_times = sol_general.t
# site_names = distance_df.site
site_names = ["Porto Moniz", "Pacl do Mar", "Funchal", "Desertas", "Canidal", "Santa Cruz", "Ribeira Brava", "So Vicente"]
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


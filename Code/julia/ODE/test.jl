using Pkg
Pkg.activate(".")
using OrdinaryDiffEq
using CairoMakie


# u0 = [[1_800.0 for j in 1:2] for i in 1:8]
# u0 = reduce(vcat, u0)
u0 = reshape([1_800.0 for j in 1:24], 8, 3)
r = [10.0, 0.99, 0.92]
d = [0.001, 0.001, 0.001]
K = 40_000
α = [0.1, 0.1, 0.1]
p = [r, d, K, α]

function t1!(du, u, p, t)
  counter = 0
  for site in 1:8
    counter += 1
    stage = 1
    prev_stage = 2
    du[counter] = (p[1][stage] * u[site, prev_stage] * ((p[3] - u[site, prev_stage]) / p[3])) -
                  (p[1][stage+1] * u[site, stage]) -
                  (p[2][stage] * u[site, stage])
  end
end


## defining model parameters
function exploit(t, rate)
  if (t % 365) / 365 < 0.42
    return rate
  else
    return 0.0
  end
end

function reproductive_cycle(t)
  if (t % 365) / 365 >= 0.42
    return 1.0
  else
    return 0.0
  end
end

function t2!(du, u, p, t)
  nsites = 8
  nstages = 3

  counter = 0
  for site in 1:nsites

    #stage 1
    counter += 1
    stage = 1
    prev_stage = nstages
    dispersal_probs1 = fill(0.1, nsites)
    dispersal_probs2 = fill(0.09, nsites)
    du[counter] = (reproductive_cycle(t) * p[1][stage] * u[site, prev_stage] * ((p[3] - u[site, prev_stage]) / p[3])) -
                  (p[1][stage+1] * u[site, stage]) -
                  (p[2][stage] * u[site, stage]) +
                  (sum(dispersal_probs1 .* u[:, stage])) -
                  (u[site, stage] * sum(dispersal_probs2))

    counter += 1
    stage = 2
    prev_stage = 1
    dispersal_probs1 = fill(0.1, nsites)
    dispersal_probs2 = fill(0.09, nsites)
    du[counter] = (p[1][stage] * u[site, prev_stage]) -
                  (p[1][stage+1] * u[site, stage]) -
                  (p[2][stage] * u[site, stage]) +
                  (sum(dispersal_probs1 .* u[:, stage])) -
                  (u[site, stage] * sum(dispersal_probs2))

    counter += 1
    stage = 3
    prev_stage = 2
    dispersal_probs1 = fill(0.1, nsites)
    dispersal_probs2 = fill(0.09, nsites)
    du[counter] = (p[1][stage] * u[site, prev_stage]) -
                  (p[1][stage] * u[site, stage]) -
                  (p[2][stage] * u[site, stage]) +
                  (sum(dispersal_probs1 .* u[:, stage])) -
                  (u[site, stage] * sum(dispersal_probs2))

  end
end

tspan = (100.0, 200.0)
prob = ODEProblem(t2!, u0, tspan, p)
sol = solve(prob, Rosenbrock23(), dt=0.0001, adaptive=false)

fig = Figure()
ax1 = Axis(fig[1, 1])
stage = 1
lines!(ax1, sol.t, [sol.u[i][1, stage] for i in 1:length(sol.t)], yscale=:log10)
for site in 2:8
  lines!(ax1, sol.t, [sol.u[i][site, stage] for i in 1:length(sol.t)])
end
save("figs/test_stage=$stage.pdf", fig)


### Start with a simple model
################################################################

using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using CairoMakie
using Statistics
using DataFrames
using CSV


## Starting the model
#---------------------------
nsites = 1
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

"""
captures the reproductive cycle and equals to false when fishing and to true when organisms reproduce
"""
function reproductive_cycle(t)
  if (t % 365) / 365 >= 0.42
    return true
  else
    return false
  end
end

# conversion rates between stages
# average oocytes per year per adult
avg_oocytes = mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
reggs = reggs / 500 # because the rate is too high to be handled
r = [reggs, 0.998611, 0.971057, 0.683772, 0.00629]
# natural death rates per life stage.
d = [0.001, 0.001, 0.001, 0.001, 0.000322]
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


function nofishing!(du, u, p, t)
  # Change these parameters if you change the model
  nsites = 1
  nlifestages = 5

  r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α = p

  counter = 0
  for site in 1:nsites

    # Stages 1 Egg
    counter += 1
    stage = 1
    prev_stage = 5

    du[counter] = (r[stage] * u[site, prev_stage]) -
                  (r[stage+1] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    # Stage 2 Trochophore
    counter += 1
    stage = 2
    prev_stage = 1
    du[counter] = (r[stage] * u[site, prev_stage]) -
                  (r[stage+1] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    #stage 3 Veliger
    counter += 1
    stage = 3
    prev_stage = 2
    du[counter] = (r[stage] * u[site, prev_stage] * ((K - u[site, stage]) / K)) -
                  (r[stage+1] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    # stage 4 Juvenile
    counter += 1
    stage = 4
    prev_stage = 3
    du[counter] = (r[stage] * u[site, prev_stage]) -
                  (r[stage+1] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    # stage 5 adult
    counter += 1
    stage = 5
    prev_stage = 4
    du[counter] = (r[stage] * u[site, prev_stage]) -
                  (0.0 * u[site, stage]) -
                  (d[stage] * u[site, stage])

    # adult sizes
    # adult sizes. dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
    counter += 1
    du[counter] = size_growth_rate * u[site, nlifestages+1] * (1 - u[site, nlifestages+1] / (size_max - (size_max * 0.0)))

  end
end

function fishing!(du, u, p, t)
  # Change these parameters if you change the model
  nsites = 1
  nlifestages = 5
  
  r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α = p
  counter = 0
  for site in 1:nsites

    # Stages 1 Egg
    counter += 1
    stage = 1
    prev_stage = 5

    du[counter] = (0.0 * r[stage] * u[site, prev_stage]) -
                  (r[stage+1] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    # Stage 2 Trochophore
    counter += 1
    stage = 2
    prev_stage = 1

    du[counter] = (r[stage] * u[site, prev_stage]) -
                  (r[stage+1] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    #stage 3 Veliger
    counter += 1
    stage = 3
    prev_stage = 2

    du[counter] = (r[stage] * u[site, prev_stage] * ((K - u[site, stage]) / K)) -
                  (r[stage+1] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    # stage 4 Juvenile
    counter += 1
    stage = 4
    prev_stage = 3

    du[counter] = (r[stage] * u[site, prev_stage]) -
                  (r[stage+1] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    # stage 5 adult
    counter += 1
    stage = 5
    prev_stage = 4
    du[counter] = (r[stage] * u[site, prev_stage]) -
                  (exploitation_rates[site] * u[site, stage]) -
                  (d[stage] * u[site, stage])

    # adult sizes
    # adult sizes. dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
    counter += 1
    du[counter] = size_growth_rate * u[site, nlifestages+1] * (1 - u[site, nlifestages+1] / (size_max - (size_max * exploitation_rates[site])))

  end
end

function split_points(t)
  times = t[1]:0.01:t[2]
  values = [reproductive_cycle(i) for i in times]
  periods = [[[times[1]], values[1]]]
  for v in 2:length(values)
    if values[v] != periods[end][2]
      push!(periods[end][1], times[v-1])
      push!(periods, [[times[v]], values[v]])
    end
  end
  push!(periods[end][1], t[2])
  return periods
end

function combined(u, p, t)
  periods = split_points(t)
  all_u = Matrix{Float64}[]
  all_times = Float64[]
  for period in periods
    func = period[2] ? nofishing! : fishing!
    prob = ODEProblem(func, u, period[1], p)
    sol = solve(prob)
    all_u = vcat(all_u, sol.u)
    all_times = vcat(all_times, sol.t)
    u = sol[end]
  end

  return all_u, all_times
end

tspan_general = (1.0, 3000.0)
all_u, all_times = combined(u0_general, p_general, tspan_general)


# Plot
# site_names = distance_df.site
site_names = ["Porto Moniz"]#, "Pacl do Mar", "Funchal"]#, "Desertas", "Canidal", "Santa Cruz", "Ribeira Brava", "So Vicente"]
site_indices = [0]#, 6, 12]#, 18, 24, 30, 36, 42]#, 48, 54]
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
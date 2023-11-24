using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using JLD2
include("load_params.jl")
exploitation_rates = [0.9880685477383748, 0.9880947574980841, 0.9891331673866427, 0.9874718910553065, 0.9875900832783625, 0.9900632400512024, 0.9894700088322096, 0.9894645471975891]
exploitation_rates_org = exploitation_rates .- exploitation_rates[4]

e_factors = [0.1, 0.5, 1.0, 2.0, 10.0, 50.0, 80.0, 100.0]

for factor in e_factors
  exploitation_rates = exploitation_rates_org * factor
  p_general = [r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α, mig_probs]
  p_general_flat = flatten_p(p_general)

  start_year = 2006
  end_year = 2018
  nyears = end_year - start_year + 1
  ndays = nyears * 365.0
  tspan_general = (0.0, ndays)

  prob_general = ODEProblem(nsites!, u0_general, tspan_general, p_general_flat)
  sol_general = solve(prob_general)

  # Saving the results
  df_general = DataFrame(sol_general)
  CSV.write("../figs/general_e_factor=$(factor).csv", df_general)

  # Plotting the results
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
    save("../figs/stage=$(stage)_e_factor=$(factor).pdf", fig)
  end

  # Changes in body size
  stage = 6
  fig, ax, plt = lines(all_times, [all_u[i][1, stage] for i in 1:length(all_times)], yscale=:log10, label=site_names[1])
  ax.title = "Body size"
  for site in 2:nsites
    lines!(ax, all_times, [all_u[i][site, stage] for i in 1:length(all_times)], label=site_names[site])
  end
  fig[1, 2] = Legend(fig, ax, "Site")
  save("../figs/body_sizes_e_factor=$(factor).pdf", fig)
end

########################################################
## Plot mean population size vs exploitation rate
########################################################

site_names = ["Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente"]
lifestage_names = ["Eggs", "Trochophore", "Veliger", "Juvenile", "Adult", "Adult Body Size"]

fig_all = Figure(resolution=(600, 2400))

for lifestage in 1:6
  # Variable to hold each site's mean population size for each e_factor
  sites_mean_pop_size = zeros(Float64, 8, length(e_factors))

  # Iterate through e_factors
  for (i, e_factor) in enumerate(e_factors)
      filename = "../figs/general_e_factor=$(e_factor).csv"

      # Load the CSV file into a DataFrame
      df = CSV.File(filename) |> DataFrame

      # Iterate through sites
      for site in 1:8
          # Create index for the site's relevant column in the flat matrix
          index = ((lifestage-1) * 8 + site) + 1  # +1 to account for the time column

          # Calculate the mean population size for the site and store it
          sites_mean_pop_size[site, i] = mean(df[:, index])
      end
  end

  # Save the plot for each lifestage separately.
  # Create the line plot
  fig = Figure(resolution = (600, 400))
  ax = Axis(fig[1, 1], xlabel = "Exploitation rate", ylabel = "Mean Population Size", title = lifestage_names[lifestage])

  for site in 1:8
    lines!(ax, e_factors, sites_mean_pop_size[site, :], label=site_names[site])
  end

  # Add a legend
  leg = Legend(fig[1, 2], ax)
  fig[1, 1] = ax
  # axislegend(ax, "Site $(1:8)", position = :rt)

  save("../figs/exploitation_vs_population_size_lifestage=$(lifestage).pdf", fig)

  # Save all life stages in the same plot
  # Create the line plot in a new subplot
  ax = Axis(fig_all[lifestage, 1], xlabel="Exploitation rate", ylabel="Mean Population Size", title=lifestage_names[lifestage])

  for site in 1:8
    lines!(ax, e_factors, sites_mean_pop_size[site, :], label=site_names[site])
  end

  # Add a legend
  leg = Legend(fig_all[lifestage, 2], ax)
  fig_all[lifestage, 1] = ax
end

# Save the figure
save("../figs/exploitation_vs_population_size_all_lifestages.pdf", fig_all)

########################################################
## Sensitivity analysis
########################################################

# for lifestage in 1:6
lifestage=5

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
isfile(output) ? m = load(output, "results") : m = gsa(fn, Sobol(), bounds_n, samples=100) && save(output, "results", m)

# Plot the results
f = Figure(resolution=(600, 2400))
for site in 1:8
  param_names = ["r", "d", "size_growth_rate", "distance_matrix", "exploitation_rates", "size_max", "K", "α", "mig_probs"]
  important_parameter_indices = vcat(collect(1:5), # r for all life stages
    vcat(6:10), # d for all life stages
    vcat(76:83), # exploitation rates for all sites. (10+64+1):(10+64+1+8)
    vcat(86:87), # α
  )
  important_parameter_names = ["r1", "r2", "r3", "r4", "r5", "d1", "d2", "d3", "d4", "d5", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "α1", "α2"]  ## Bar plot
  normalized_mST = m.ST ./ sum(m.ST, dims=2)
  sensitivity_indices = normalized_mST[site, important_parameter_indices]
  if site == 4
    ax  = Axis(f[site, 1], xlabel="Parameter", ylabel="Sensitivity index - total effect", xticks=(1:length(important_parameter_names), important_parameter_names), title="$(site_names[site])")
  else
    ax  = Axis(f[site, 1], xlabel="Parameter", ylabel="Sensitivity index - total effect", xticks=(1:length(important_parameter_names), important_parameter_names), title="$(site_names[site])", yscale = Makie.pseudolog10)
  end
  ylims!(ax, 0.0, 1.0)

  bars = barplot!(ax, 1:length(important_parameter_names), sensitivity_indices, color=:black)
  # site != 8 &&  hidexdecorations!(ax, grid=true, ticks=true)
end
save("../figs/sensitivity_indices_barplot_total_effect_lifestage=$(lifestage)_e_factor=$(factor).pdf", f)
# end
# end
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
exploitation_rates_org = cat(exploitation_rates, exploitation_rates, dims=2) # for two species

e_factors = [0.1, 0.5, 1.0, 2.0, 10.0, 50.0, 80.0, 100.0]
competition_coefficients = [[0.5, 0.5], [1.0, 1.0]]
species = ["P. ordinaria", "P. aspera"]
nspecies = length(species)

for cf in competition_coefficients, ef in e_factors
  exploitation_rates = exploitation_rates_org * ef
  p_general = [r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α, mig_probs, cf]
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
  CSV.write("../figs/general_e_factor=$(ef)_competition_coef=$(cf).csv", df_general)
  @save "../figs/general_e_factor=$(ef)_competition_coef=$(cf).jld2" sol_general

  # Plotting the results
  ## Load the ODE results
  # JLD2.@load "../figs/general.jld2" sol_general

  all_u = sol_general.u
  all_times = sol_general.t
  # site_names = distance_df.site
  site_names = ["Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente"]
  colors = ["#55a47b", "#a361c7", "#91db6f", "#e87cb2", "#a6953f", "#6588cd", "#ffa66a", "#cb5358"]
  linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]  # Define more styles if you have more species

  for stage in 1:5
    fig = Figure(resolution=(900, 600), backgroundcolor=:transparent)
    ax1 = Axis(fig[1, 1], backgroundcolor=:transparent)

    legends_colors = []  # Store legends for each site
    legends_styles = []  # Store legends for each species
    plots_colors = []  # Store plots for color legend
    labels_colors = []  # Store labels for color legend
    plots_styles = []  # Store plots for style legend
    labels_styles = []  # Store labels for style legend

    for sp in 1:nspecies
      for site in 1:nsites
        p = lines!(ax1, all_times, [all_u[i][site, stage, sp] for i in 1:length(all_times)], yscale=:log10, color=colors[site], linestyle=linestyles[sp])
        ax1.title = "N for stage: $(stage)"
        ax1.xlabel = "Time (days)"
        ax1.ylabel = "Population size"
        if sp == 1  # Only add to color legend for the first species
          push!(plots_colors, p)
          push!(labels_colors, site_names[site])
        end
      end
      # Create a legend for each species with a dummy plot for the line style
      p_dummy = lines!(ax1, [NaN], [NaN], linestyle=linestyles[sp])  # Dummy plot, won't be visible
      push!(plots_styles, p_dummy)
      push!(labels_styles, species[sp])
    end

    # Create legends
    push!(legends_colors, Legend(fig, plots_colors, labels_colors, "Site", labelsize=20))
    push!(legends_styles, Legend(fig, plots_styles, labels_styles, "Species", labelsize=20))

    fontsize_theme = Theme(fontsize=20)
    set_theme!(fontsize_theme)

    # Combine the legends
    fig[1, 2] = Makie.vgrid!(legends_colors...)
    fig[1, 3] = Makie.vgrid!(legends_styles...)

    save("../figs/stage=$(stage)_e_factor=$(ef)_competition_coef=$(cf).pdf", fig)
    save("../figs/stage=$(stage)_e_factor=$(ef)_competition_coef=$(cf).png", fig)

  end

  # Changes in body size
  stage = 6

  fig = Figure(resolution=(900, 600), backgroundcolor=:transparent)
  ax2 = Axis(fig[1, 1], backgroundcolor=:transparent)

  legends_colors = []  # Store legends for each site
  legends_styles = []  # Store legends for each species
  plots_colors = []  # Store plots for color legend
  labels_colors = []  # Store labels for color legend
  plots_styles = []  # Store plots for style legend
  labels_styles = []  # Store labels for style legend

  for sp in 1:nspecies
    for site in 1:nsites
      p = lines!(ax2, all_times, [all_u[i][site, stage, sp] for i in 1:length(all_times)], yscale=:log10, color=colors[site], linestyle=linestyles[sp])
      ax2.title = "Body size"
      ax2.xlabel = "Time (days)"
      ax2.ylabel = "Body size (mm)"
      if sp == 1  # Only add to color legend for the first species
        push!(plots_colors, p)
        push!(labels_colors, site_names[site])
      end
    end
    # Create a legend for each species with a dummy plot for the line style
    p_dummy = lines!(ax2, [NaN], [NaN], linestyle=linestyles[sp])  # Dummy plot, won't be visible
    push!(plots_styles, p_dummy)
    push!(labels_styles, species[sp])
  end

  # Create legends
  push!(legends_colors, Legend(fig, plots_colors, labels_colors, "Site", labelsize=20))
  push!(legends_styles, Legend(fig, plots_styles, labels_styles, "Species", labelsize=20))

  fontsize_theme = Theme(fontsize=20)
  set_theme!(fontsize_theme)

  # Combine the legends
  fig[1, 2] = Makie.vgrid!(legends_colors...)
  fig[1, 3] = Makie.vgrid!(legends_styles...)

  save("../figs/body_size_e_factor=$(ef)_competition_coef=$(cf).pdf", fig)
  save("../figs/body_size_e_factor=$(ef)_competition_coef=$(cf).png", fig)

end

########################################################
## Plot mean population size vs exploitation rate TODO
########################################################

site_names = ["Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente"]
lifestage_names = ["Eggs", "Trochophore", "Veliger", "Juvenile", "Adult", "Adult Body Size"]

fig_all = Figure(resolution=(600, 2400))

for lifestage in 1:6
  # Variable to hold each site's mean population size for each e_factor
  sites_mean_pop_size = zeros(Float64, 8, length(e_factors), nspecies)

  # Iterate through e_factors
  for (i, ef) in enumerate(e_factors), (j, cf) in enumerate(competition_coefficients)
      filename = "../figs/general_e_factor=$(ef)_competition_coef=$(cf).csv"

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
  colors = ["#55a47b", "#a361c7", "#91db6f", "#e87cb2", "#a6953f", "#6588cd", "#ffa66a", "#cb5358"]

  fig = Figure(resolution=(900, 600), backgroundcolor=:transparent)
  ax = Axis(fig[1, 1], xlabel = "Exploitation rate", ylabel = "Mean Population Size", title = lifestage_names[lifestage], backgroundcolor=:transparent)

  for site in 1:8
    lines!(ax, e_factors, sites_mean_pop_size[site, :], label=site_names[site], color=colors[site])
  end
  fontsize_theme = Theme(fontsize=20)
  set_theme!(fontsize_theme)
  # Add a legend
  leg = Legend(fig[1, 2], ax, labelsize=20)
  fig[1, 1] = ax
  # axislegend(ax, "Site $(1:8)", position = :rt)

  save("../figs/exploitation_vs_population_size_lifestage=$(lifestage).pdf", fig)
  save("../figs/exploitation_vs_population_size_lifestage=$(lifestage).png", fig)

  # Save all life stages in the same plot
  # Create the line plot in a new subplot
  ax = Axis(fig_all[lifestage, 1], xlabel="Exploitation rate", ylabel="Mean Population Size", title=lifestage_names[lifestage], backgroundcolor=:transparent)

  for site in 1:8
    lines!(ax, e_factors, sites_mean_pop_size[site, :], label=site_names[site], color=colors[site])
  end
  fontsize_theme = Theme(fontsize=20)
  set_theme!(fontsize_theme)
  # Add a legend
  leg = Legend(fig_all[lifestage, 2], ax, labelsize=20)
  fig_all[lifestage, 1] = ax
end

# Save the figure
save("../figs/exploitation_vs_population_size_all_lifestages.pdf", fig_all)
save("../figs/exploitation_vs_population_size_all_lifestages.png", fig_all)


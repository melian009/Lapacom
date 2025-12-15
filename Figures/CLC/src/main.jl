using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
# using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using JLD2
# using DiffEqParamEstim
# using Optim

include("load_params.jl")

start_year = 2006
end_year = 2018
nyears = end_year - start_year + 1
ndays = nyears * 365.0
tspan_general = (0.0, ndays)

prob_general = ODEProblem(nsites!, u0_general, tspan_general, p_general_flat)
sol_general = solve(prob_general)

# Saving the results
@save "../figs/general.jld2" sol_general

# Plotting the results
## Load the ODE results
JLD2.@load "../figs/general.jld2" sol_general

all_u = sol_general.u
all_times = sol_general.t
# site_names = distance_df.site
site_names = ["Porto Moniz", "Paúl do Mar", "Funchal", "Desertas", "Caniçal", "Santa Cruz", "Ribeira Brava", "São Vicente"]
# colors = cgrad(:lightrainbow, length(site_names), categorical=true)
colors = ["#55a47b", "#a361c7", "#91db6f", "#e87cb2", "#a6953f", "#6588cd", "#ffa66a", "#cb5358"]

for stage in 1:5
  fig = Figure(resolution=(900, 600), backgroundcolor=:transparent)
  ax1 = Axis(fig[1, 1], backgroundcolor=:transparent)
  lines!(ax1, all_times, [all_u[i][1, stage] for i in 1:length(all_times)], yscale=:log10, label=site_names[1], color=colors[1])
  ax1.title = "N for stage: $(stage)"
  for site in 2:nsites#2:8
    lines!(ax1, all_times, [all_u[i][site, stage] for i in 1:length(all_times)], label=site_names[site], color=colors[site])
  end
  fontsize_theme = Theme(fontsize=20)
  set_theme!(fontsize_theme)
  fig[1, 2] = Legend(fig, ax1, "Site", labelsize=20)
  save("../figs/stage=$stage.pdf", fig)
  save("../figs/stage=$stage.png", fig)
end

# Changes in body size
stage = 6
fig = Figure(resolution=(900, 600), backgroundcolor=:transparent)
ax1 = Axis(fig[1, 1], backgroundcolor=:transparent)
lines!(ax1, all_times, [all_u[i][1, stage] for i in 1:length(all_times)], yscale=:log10, label=site_names[1], color=colors[1])
ax1.title = "Body size"
for site in 2:nsites
  lines!(ax1, all_times, [all_u[i][site, stage] for i in 1:length(all_times)], label=site_names[site], color=colors[site])
end
fontsize_theme = Theme(fontsize=20)
set_theme!(fontsize_theme)
fig[1, 2] = Legend(fig, ax1, "Site", labelsize=20)
save("../figs/body_sizes.pdf", fig)
save("../figs/body_sizes.png", fig)
using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
# using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
# using DiffEqParamEstim
# using Optim

include("load_params.jl")

tspan_general = (0.0, 1000.0)
prob_general = ODEProblem(nsites!, u0_general, tspan_general, p_general_flat)
sol_general = solve(prob_general)


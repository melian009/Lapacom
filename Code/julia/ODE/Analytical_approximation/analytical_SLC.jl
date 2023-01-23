using Pkg
using ForwardDiff
#Pkg.activate(".")
using LinearAlgebra
# using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
#using CairoMakie
using Statistics
using DataFrames
# using CSV
using DiffEqParamEstim
using Optim
using Plots
using Plots.PlotMeasures

import ForwardDiff.jacobian


function f(x)
    F = zero.(x)
    
    
    
    function single_site_S!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = p
  du[1] = dNⱼ = (r * Nₐ * (Sₐ/sizeₘₐₓ) * ((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
end

    F[1] = x[1]^2 + x[3]
    F[2] = x[1] + x[2]
    F[3] = x[2]^2 + x[3]^2

    return F
end

x0 = [1,2,3]
J0 = jacobian(f, x0)

#Output
#3×3 Matrix{Int64}:
# 2  0  1
# 1  1  0
# 0  4  6


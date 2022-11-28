## Packeges
using Pkg
Pkg.activate(".")
using LinearAlgebra
#using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Plots

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
  du[1] = dNⱼ = (r * Nₐ * (Sₐ/sizeₘₐₓ)*((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
end

Et(t) = (sin(t)^2) / 2  # time varying exploitation
# To adjust the sin function horizonally (stretch/shrink the wave length), multiply t by a factor c.
# You may also use a function with if statements.
function Et(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return 0.5
  end
end

p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0]
u0_1 = [1e3, 1e3, 40.0]
tspan_1 = (0.0, 200.0)
prob_1 = ODEProblem(single_site!, u0_1, tspan_1, p_1)
sol_1 = solve(prob_1, Tsit5())

plot(sol_1,vars=(0, 1), label="Nⱼ")
plot!(sol_1,vars=(0, 2), label="Nₐ")
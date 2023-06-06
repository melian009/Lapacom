#Packages
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
using Symbolics
import ForwardDiff.jacobian

function SLC!(du, u, p, t)
   Na, Sa = u
   i, r, K, H, X, da, Smax, gamma = p
   Saverage = du[2]
   Smaturity = calculate_size_at_first_maturity(Saverage)

  du[1] = dNa = X(t) * r * Na * Rep_cap(Saverage, Smaturity, Smax) * (K - Na/K) - (1 - X(t)) * H(i) * Na - (da[i] * Na) 
  du[2] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 - H[i])))
end



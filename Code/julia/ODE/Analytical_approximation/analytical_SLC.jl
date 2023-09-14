using Pkg
Pkg.activate(".")
using ForwardDiff
using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
using Statistics
using DiffEqParamEstim
using Optim
using Plots
using Plots.PlotMeasures
using Symbolics
import ForwardDiff.jacobian
using NonlinearSolve


@variables Na Sa r X R d K H g Smax

J_SLC = Symbolics.jacobian([(X * r * R * Na)*((K - Na)/K) - (d * Na) - H*(1 - X)*Na, (g * Sa) * (1 - (Sa/(Smax * (1 - H * (1-X)))))], [Na, Sa])

#OUTPUT
# [Na, Sa])
#2×2 Matrix{Num}:
# (R*X*r*(K - Na) - Na*R*X*r) / K - d - H*(1 - X)                                                                          0
#                                               0  g*(1 + (-Sa) / (Smax*(1 - H*(1 - X)))) + (-Sa*g) / (Smax*(1 - H*(1 - X)))

det_jac = det(J_SLC)

#OUTPUT
#(g*(1 + (-Sa) / (Smax*(1 - H*(1 - X)))) + (-Sa*g) / (Smax*(1 - H*(1 - X))))*((R*X*r*(K - Na) - Na*R*X*r) / K - d - H*(1 - X))
   
M = Symbolics.simplify(det_jac)

#OUTPUT
#((Smax*g + H*Smax*X*g - 2.0Sa*g - H*Smax*g)*(H*K*X + K*R*X*r - H*K - K*d - 2Na*R*X*r)) / (K*Smax*(1 + H*X - H))

Symbolics.solve_for([(X * r * R * Na)*((K - Na)/K) - d * Na - H*(1 - X)*Na, g * Sa * (1 - (Sa/(Smax * (1 - H * (1-X)))))], [Na,Sa]; simplify=true)

#OUTPUT

#EXAMPLE
#https://github.com/JuliaSymbolics/Symbolics.jl/issues/381
#@variables r y x₂ h₁ g₁ x₁ g₂ h₂ g₄ g₃ f;
#Symbolics.solve_for([(r - y - x₂ * h₁) * g₁ ~ x₁,(x₁ - h₂ * y) * g₂ ~ x₂,x₂ * g₃ + x₁ * g₄ + f ~ y], [x₁,x₂,y]; simplify=false)


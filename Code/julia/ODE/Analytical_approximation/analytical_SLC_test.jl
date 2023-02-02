using Pkg
using ForwardDiff
#Pkg.activate(".")
using LinearAlgebra
# using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
#using CairoMakie
using Statistics
#using DataFrames
# using CSV
using DiffEqParamEstim
using Optim
using Plots
using Plots.PlotMeasures

import ForwardDiff.jacobian


#Test 1

function f(x)
    F = zero.(x)
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

#Test 2
Symbolics.jacobian([x + x*y, x^2 + y],[x, y])




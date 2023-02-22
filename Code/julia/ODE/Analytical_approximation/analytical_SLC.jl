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
using Symbolics
using SymPy
import ForwardDiff.jacobian


@variables x[1] x[2]
@parameters r K g d[1] d[2] E

function f(x)
    F = zero.(x)
    F[1]  = (r * x[2] * ((K - x[2]) / K)) - (d[1] * x[1]) - (g * x[1])
    F[2]  = (g * x[1]) - (d[2] * x[2]) - (E * x[2])
    return F
end

x0 = [1,2]
J0 = jacobian(f, x0)

#Output
#3×3 Matrix{Int64}:
# 2  0  1
# 1  1  0
# 0  4  6




#= Set up symbolic variables and parameters
@variables x[1] x[2]
@parameters r K g d[1] d[2] E

# D = Differential(t)
F[1]  = (r * x[2] * ((K - x[2]) / K)) - (d[1] * x[1]) - (g * x[1])
F[2]  = (g * x[1]) - (d[2] * x[2]) - (E * x[2])


  ## Inputs
# Make Lagrangian
V = -mₚ*Lₚ*g*cos(θ)
T = (mₖ+mₗ)*v^2/2 + mₚ*(Lₚ^2*ω^2 + 2*Lₚ*ω*v*cos(θ) +v^2)/2
L = T - V

# Cart input force
F = 1000sin(t)

# Generalized forces
Q = [F, 0]

# Make equations of motion
slosh_cart = LagrangeEOM(L, [v, ω], [x, θ], [Lₐ, Lₚ, mₖ, mₗ, mₚ, g], t; Q)

# Initial Conditions
ic = [
    θ => deg2rad(10)
]

# Parameters
p = Dict([
    Lₐ => 10
    Lₚ => 0.5
    mₖ => 100
    mₗ => 0
    mₚ => 25
    g => 9.80665
])


## Simulation
sol = solve(ODEProblem(slosh_cart, ic, (0.0, 10.0), [p...]))
=#



@variables x1 x2 r K g d1 d2 E
#Symbolics.jacobian([y2 + y2*y1, y1^2 + y1],[y1, y2])
J = Symbolics.jacobian([(r * x2 * ((K - x2) / K)) - (d1 * x1) - (g * x1),
(g * x1) - (d2 * x2) - (E * x2)],[x1, x2])
Det_J = det(J)



M = Symbolics.simplify(Det_J)



Exp_lim = 1                 # Exploitation max limit 
m=0.05                      # Interval of exploitation values 
Expl= 0:m:Exp_lim 
X_2= zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
c=0
for n = 0:m:Exp_lim
r = 0.6
g = 0.06
d1 = 0.05
d2 = 0.08 
E = n 
K=1e4 
x2= K*(g*(r-E-d2) - d1*(E + d2))/(2*g*r)
c=c+1
X_2[c,] = x2
end

plot(Expl,X_2)

#= This is the inecuality for the adult abundances
 Equation:
 X2= K * ( g * (r - E - d2) - d1 * (E + d2) ) / (2 * g * r)
 Conditions:
 X2>0 if  g * (r - E - d2)  > d1 * (E + d2)
=#

#Packages
using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using GlobalSensitivity
#using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Symbolics
#using GLMakie
using Plots
import ForwardDiff.jacobian

#Check SLC.jl 

function my_ode!(du, u, t, p)
    x1, x2, y1, y2 = u
    re1, re2, K, H1, H2, X, R1, R2, g1, g2, d1, d2, c12, c21, Smax = p
    
    #Parameter values species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)
    # Population Growth rate estimation (r=reggs):

   # oocytes_po = 385613                  # Average: Patella ordinaria (nº of Eggs)
   # oocytes_pa = 73029                   # Average: Patella aspera (nº of Eggs)
   # oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
   # reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.
   # re = reggs / 500     # Population growth rate: re[1,1] re[2,1]
    
   R1 = 1
   R2 = 1
   X =  0
   K = 640000          # Carrying capacity 
   H1 = 0.639
   H2 = 0.57 # Exploitation rate (H)
   g1 = 0.32
   g2 = 0.36 # Adult growth rate
    #rates2 = [0.02,0.01]
    #gEA = 0.006          # Instant conversion between stages.
   d1 = 0.55
   d2 = 0.59    # Natural mortality rate for adults
   Smax = 56              # Maximum size for adults
   re1 = 0.32
   re2 = 0.36 # Adult growth rate
   c12 = 0.05#competition term species 2 on 1
   c21 = 0.05#competition term species 1 on 2
    
    du[1] = (X * re1 * R1 * x1)*((K - x1)/K) - (d1 * x1) - H1*(1 - X)*x1 - c12*x1*x2
    du[2] = (X * re2 * R2 * x2)*((K - x2)/K) - (d2 * x2) - H2*(1 - X)*x2 - c21*x1*x2
    du[3] = (g1 * y1) * (1 - (y1/(Smax * (1 - H1 * (1-X)))))
    du[4] = (g2 * y2) * (1 - (y2/(Smax * (1 - H2 * (1-X)))))
end

u0 = [100.0, 100.0, 25.0, 25.0]  # Initial conditions
tspan = (0.0, 1000.0)  # Time span for the simulation (from t=0 to t=1000)
solver = Tsit5()
prob = ODEProblem(my_ode!, u0, tspan)
sol = solve(prob, solver) # Eros: There is an error when I run this line and I don't know how to fix it. "Error: BoundsError: attempt to access Float64 at index [2]"

plot(sol, xlabel="Time", ylabel="State Variables", label=["x1" "x2" "y1" "y2"])




using Pkg
Pkg.activate(".")
using ForwardDiff
using LinearAlgebra
# using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
# using CairoMakie
using Statistics
# using DataFrames
# using CSV
using DiffEqParamEstim
using Optim
using Plots
using Plots.PlotMeasures
using Symbolics
import ForwardDiff.jacobian

#=
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

=#


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

#======================================================================
Analytical approach for the simple life cycle in a single site (SLC-OS)
======================================================================#

```
Simple life cycle equations:
 dNe/dt = (X * r * Na *(Sa/Smax)*((K - Na)/K)) - (de * Ne) - (g * Ne)
 dNa/dt = (g * Ne) - (da * Na) - (Na * H * (1-X))
 dSa/dt = size_growth_rate * Sa * (1 - Sa/(Smax * (1 - H * (1-X))))
```

@variables Na Ne Sa r K de da g X H Smax size_growth_rate 

# Symbolics.jacobian([f1(y1,y2), f2(y1,y2)],[y1, y2])

J_SLC = Symbolics.jacobian([
 (X * r * Na *(Sa/Smax)*((K - Na)/K)) - (de * Ne) - (g * Ne),
 (g * Ne) - (da * Na) - (Na * H * (1-X)),
 size_growth_rate * Sa * (1 - Sa/(Smax * (1 - H * (1-X))))], # = dSa/dt
[Ne, Na, Sa]) #Vairables a considerar para calcular la matriz jacobiana  


#J_SLC = [(-de-g) ((Sa*r*(K-2 * Na))/(K*Smax)) ((Na*r*(K-Na))/(K*Smax));g (-E-da) 0;0 0 (size_growth_rate*Sa*(1-(2*Na)/(Smax*(1-E))))]

# Cálculo del determinante
Det_SLC = det(J_SLC)



#Simplificación del determinante para despejar las las variables.
M = Symbolics.simplify(Det_SLC)


# Symbolics.solve_for(M, Sa) # Error

Na = Symbolics.solve_for(M, Na) #

dNedt = (X * r * Na *(Sa/Smax)*((K - Na)/K)) - (de * Ne) - (g * Ne)
dNadt = (g * Ne) - (da * Na) - (Na * H * (1-X))
dSadt = size_growth_rate * Sa * (1 - Sa/(Smax * (1 - H * (1-X))))


Ne = Symbolics.solve_for(dNadt, Ne)

Sa = Symbolics.solve_for(dNadt, Sa)

#Estimamos los valores de Na y Sa para distintos valores de E.

Exp_lim = 1                 # Exploitation max limit 
m=0.005                      # Interval of exploitation values 
Expl= 0:m:Exp_lim           # Exploitation ranges
Sa_ = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
Na_ = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
c=0

for n = 0:m:Exp_lim
  function Et(t)
  if modf(t)[1] < 0.5         # 50% of the adult population is exploited
    return 0.0
  else
    return n                    #For an exploitation value equal to 1, the mathematical result is erratic because the size equation will present a denominator division equal to 0
  end
 end
 p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40] # r, g, dⱼ, dₐ, E, K, size_growth_rate, Smax
 u0 = [1e3,40]
 tspan = [0,365*2]
 prob_1 = ODEProblem(analitical_aproach_SLC_SS!, u0, tspan, p_1)
 sol_1 = solve(prob_1, Tsit5())
 c=c+1

 Na_[c,] = sol_1[1,end]
 Sa_[c,] = sol_1[2,end]

end

#Adult size plot by exploitation
plot(Expl,Sa_,label="NA: N (adults)")
xlims!(0,1)
xlabel!("Exploitation rate")
ylabel!("N (nº individuals)") 



````````````````````````````````````````````````````

````````````````````````````````````````````````````

#= Possible analitical solution for SLC in a single site
  Na = K*(2*E + g + da)*(3*E + g + 2*da)^-1
  Ne = (Na * (- E - g) * (K - Na)) * (g * (K - 2*Na))^-1
  Sa = (- E - g) * (- da - g) * K * Smax * (g*r*(K - 2*Na))^-1
=#

#Estimamos los valores de Na y Ne para distintos valores de E.

Exp_lim = 1                 # Exploitation max limit 
m=0.05                      # Interval of exploitation values 
Expl= 0:m:Exp_lim           # Exploitation ranges
Ne_ = zeros(Float64, size(Expl)) # Void vector to array number of eggs for diferent exploitation values
Na_ = zeros(Float64, size(Expl)) # Void vector to array number of adults for diferent exploitation values
Sa_ = zeros(Float64, size(Expl))
c=0

for n = 0:m:Exp_lim
r = 0.6
g = 0.06
de = 0.05
da = 0.08 
E = n 
K = 1e4 
#Possible Analytical solution
Na = K * (2 * E + g + da) * (3 * E + g + 2 * da)^-1
Ne = (Na * (- E - g) * (K - Na)) * (g * (K - 2 * Na))^-1
Sa = (- E - g) * (- da - g) * K * Smax * (g * r * (K - 2 * Na))^-1

c=c+1
Ne_[c,] = Ne
Na_[c,] = Na
Sa_[c,] = Sa
end

plot(Expl,Ne_,label="AA: N (eggs)")
plot!(Expl,Na_,label="AA: N (adults)")
xlims!(0,1)
ylims!(-1000,15000)
xlabel!("Exploitation rate")
ylabel!("N (nº individuals)") 

plot(Expl,Sa_,label="AA: Size (adults)")
xlims!(0,1)
ylims!(-1000,15000)
xlabel!("Exploitation rate")
ylabel!("N (nº individuals)") 

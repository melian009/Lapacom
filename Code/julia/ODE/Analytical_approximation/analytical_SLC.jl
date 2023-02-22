using Pkg
using ForwardDiff
# Pkg.activate(".")
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
# using SymPy
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


# Jacobian matrix extimation for a simple life cycle (only abundance equations).

@variables Ne Na r K g de da E

# Symbolics.jacobian([f1(y1,y2), f2(y1,y2)],[y1, y2])
J = Symbolics.jacobian([(r * Na * ((K - Na) / K)) - (de * Ne) - (g * Ne), # = dNe/dt 
    (g * Ne) - (da * Na) - (E * Na)], # = dNa/dt
[Ne, Na]) #Vairables a considerar para calcular la matriz jacobiana  

# Cálculo del determinante
Det_J = det(J)

#Simplificación del determinante
M = Symbolics.simplify(Det_J)

#Determinante simplificado = (E*K*de + E*K*g + K*da*de + K*da*g + 2Na*g*r - K*g*r) / K
#= 
Si Det(J)=0
 (E*K*de + E*K*g + K*da*de + K*da*g + 2Na*g*r - K*g*r) / K = 0
 (E*K*de + E*K*g + K*da*de + K*da*g + 2Na*g*r - K*g*r) = 0
 2Na*g*r = K*g*r - (E*K*de + E*K*g + K*da*de + K*da*g)

 Despejamos Na:
 Na = (K*g*r - (E*K*de + E*K*g + K*da*de + K*da*g))/(2*g*r)
=#

#=
 Simplificamos la ecuación:
 Na = (K *(g *(r - (de + E)) - da *(de + g)))/(2*g*r) 
    Si g*(r-(de + E) > da *(de + g) | Na > 0
=#
#=
 Si dNa/dt = 0 y Na = f(K,g,r,da,de):
    dNa/dt = (g * Ne) - (da * Na) - (E * Na) = 0

 Despejamos Ne:
 Ne = ((da + E) * Na)/g
=#


#Estimamos los valores de Na y Ne para distintos valores de E.

Exp_lim = 1                 # Exploitation max limit 
m=0.05                      # Interval of exploitation values 
Expl= 0:m:Exp_lim           # Exploitation ranges
Ne_ = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
Na_ = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
c=0
for n = 0:m:Exp_lim
r = 0.6
g = 0.06
de = 0.05
da = 0.08 
E = n 
K=1e4 
Na = (K *(g *(r - (de + E)) - da *(de + g)))/(2*g*r)
Ne = ((da + E) * Na)/g
c=c+1
Ne_[c,] = Ne
Na_[c,] = Na
end

#= 
This is the inecuality for the adult abundances
 Equations:
 Na = K * ( g * (r - E - d2) - d1 * (E + d2) ) / (2 * g * r)
 Ne = ((da + E) * Na)/g

 Condition:
 if  g * (r - E - d2)  > d1 * (E + d2) then Na > 0 & Ne > 0
=#

plot(Expl,Ne_,label="N (eggs)",colour="blue")
plot!(Expl,Na_,label="N (adults)",colour="red")
xlims!(0.0,1)
ylims!(-1000,8500)
xlabel!("Exploitation rate (E)")
ylabel!("Abundance (N)")

#Punto de inflexión en E = 0.4



#=
Realizamos la aproximación numérica al mismo sistema de ecuaciones 
 y comprobamos si coinciden los resultados:
=#

# SLC in One site: 
function single_site!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = pm
  du[1] = dNⱼ = (r * Nₐ * ((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
end

Exp_lim = 1                  # Exploitation max limit 
m=0.05                           # Interval of exploitation values 
Expl= 0:m:Exp_lim                  # Expoitation values for plotting
tspan = (0.0, 365)               # Time value 
N_et_1 = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
N_at_1 = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
c = 0                              # C is the position of the vector N_et, N_at and S_at

for n = 0:m:Exp_lim
  function Et(t)
  if modf(t)[1] < 0.42          # Exploitation only ocurr in the 42% of the year
    return 0.0
  else
    return n                    # For an exploitation value equal to 1, the mathematical result is erratic because the size equation will present a denominator division equal to 0
  end
 end
 u0 = [1e3,1e3]                # Initial conditions of N_e, N_a
 p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4] # r, g, dⱼ, dₐ, E, K
 prob_1 = ODEProblem(single_site!, u0, tspan, p_1)
 sol_1 = solve(prob_1, Tsit5())
 c=c+1

 N_et_1[c,] = sol_1[1,end]
 N_at_1[c,] = sol_1[2,end]

end

#Adult size plot by exploitation
plot(Expl,N_et_1,label="NA: N (eggs)")
plot!(Expl,N_at_1,label="NA: N (adults)")
plot!(Expl,Ne_,label="NA: N (adults)")
plot!(Expl,Na_,label="NA: N (adults)")
xlims!(0,1)
ylims!(-1000,15000)
xlabel!("Exploitation rate")
ylabel!("N (nº individuals)") 

# El punto de corte en N=0 en la Explotación difiere entre la aproximación analitica y la aproximación numérica.
# Analytical approach = 0.403
# Numérical approach = 0.55


#======================================================================
Analytical approach for the simple life cycle in a single site (SLC-OS)
======================================================================#

```
Simple life cycle equations:
 dNe/dt = (r * Na *(Sa/Smax)*((K - Na)/K)) - (de * Ne) - (g * Ne)
 dNa/dt = (g * Ne) - (da * Na) - (E * Na)
 dSa/dt = size_growth_rate * Sa * (1 - Sa/(Smax * (1-E)))
```

@variables Na Ne Sa r K de da g E Smax size_growth_rate

# Symbolics.jacobian([f1(y1,y2), f2(y1,y2)],[y1, y2])
#=
J_SLC = Symbolics.jacobian([
    (r * Na *(Sa/Smax)*((K - Na)/K)) - (de * Ne) - (g * Ne), # = dNe/dt 
    (g * Ne) - (da * Na) - (E * Na), # = dNa/dt
    size_growth_rate * Sa * (1 - Sa/(Smax * (1-E)))], # = dSa/dt
[Ne, Ne, Sa]) #Vairables a considerar para calcular la matriz jacobiana  
=#

J_SLC = [(-de-g) ((Sa*r*(K-2 * Na))/(K*Smax)) ((Na*r*(K-Na))/(K*Smax));g (-E-da) 0;0 0 (size_growth_rate*Sa*(1-(2*Na)/(Smax*(1-E))))]

# Cálculo del determinante
Det_SLC = det(J_SLC)

#Simplificación del determinante para despejar las las variables.
M = Symbolics.simplify(Det_SLC)




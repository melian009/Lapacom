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
using Symbolics
# using SymPy
# import ForwardDiff.jacobian

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

Plots.plot(Expl,Ne_,label="N (eggs)",colour="blue")
Plots.plot!(Expl,Na_,label="N (adults)",colour="red")
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
  Nⱼ, Nₐ = u
  r, g, dⱼ, dₐ, E, K = p
  du[1] = (r * Nₐ * ((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
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
Plots.plot(Expl,N_et_1,label="NA: N (eggs)")
Plots.plot!(Expl,N_at_1,label="NA: N (adults)")
Plots.plot!(Expl,Ne_,label="AA: N (eggs)")
Plots.plot!(Expl,Na_,label="AA: N (adults)")
xlims!(0,1)
ylims!(-1000,15000)
xlabel!("Exploitation rate")
ylabel!("N (nº individuals)") 
savefig("SLC_AA_NA.png")
# El punto de corte en N=0 en la Explotación difiere entre la aproximación analitica y la aproximación numérica.
# Analytical approach = 0.403
# Numérical approach = 0.55
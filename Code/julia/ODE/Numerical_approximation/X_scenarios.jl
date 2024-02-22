#Packages

using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics

#=
SLC general model
dNa/dt = (R * re * Na)*((K - Na)/K) - (d * Na) - H*(1 - X(t)) * Na - c12 * Na * x2
dSa/dt = (g * Sa) * (1 - (Sa/(Smax(1 - * H*(1 - X)))))

For X=1 Species reproduces and theres no expoitation on them.
  dNa/dt = (R * re * Na)*((K - Na)/K) - (d * Na) - c12 * Na * x2
  dSa/dt = (g * Sa) * (1 - (Sa/Smax))

Non trivial solution for the scenario of X=1
    
  Nai =  K*(R*re-d-cij*Naj)/(2*R*re)
  Sai =  Snax/2

For X=0 Species reproduces and theres expoitation on them.
  dNa/dt = (R * re * Nai)*((K - Nai)/K) - (d * Na) - H * Na - c12 * Nai * Naj
  dSa/dt = (g * Sa) * (1 - (Sa/Smax))

Non trivial solutuiin for the scenario of X=0
    Nai =  K*(R*re-d-cij*Naj)/(2*R*re)
    Sai = (Smax*(1-H))/2
=#




function SLC!(du, u, p, t)
    Na, Sa = u
     X, r, K, H, d, Smax, gamma, cij,Naj_ = p
  
     # reproductive capacity
     avg_size = du[2]
     Smat = 1.34 * (avg_size) - 28.06
     R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax - Smat)), 0.0), 1.0)
  
    du[1] = dNa = r * R_ * Na * ((K - Na)/ K) - d * Na - (1 - X) * H * Na - cij * Naj_
    du[2] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - periodX(t))))
end
  
  
  
  
  
avg_oocytes = mean([92098, 804183]) # This is the actual mean. mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r_ = reggs*0.998611*0.971057*0.4820525*0.00629
# natural death rates per life stage.
d = 0.000322/365.14  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = 0.00014898749263737575

Naj_ = 92098 * 0.998611*0.971057*0.4820525*0.00629
K_ = 64000.00    # Carrying capacity
Smax_ = 53.0             # Maximum size for adults3
  #         t_0, k,  r,  K,  H ,  d, Smax,  gamma
  
X_=1
  
H1_ = [0.639,0.57] # Exploitation rate (H)

  
cij_= 0.5
  
  p_span = [X_, r_, K_, H1_, d, Smax_, size_growth_rate,cij_,Naj_] 
  n=10  #Number of years in the simulation
  tspan = (0,365.14*n)
  U0_ = [10^4, 49.25]
  prob_ = ODEProblem(SLC!, U0_, tspan, p_span)
  solve_= solve(prob_, Tsit5())
  
  t_l = length(zeros(Float64,size(1:length(solve_.u))))
  n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))
  vars = zeros(t_l,n_v)
  time2 = zeros(t_l,1)
  
  for j in 1:length(solve_.u)
    for i in 1:2
      vars[j,i] = solve_.u[j][i]
      time2[j] = solve_.t[j]
     end 
  end
  
  Na1c0 = vars[:,1]
  
  
  
  plot!(Na1c0, label="cij=1.0", color=:red, legend=:outerright, background=nothing)
  xlims!(0,100)
  ylims!(0,7*10^4)
  xlabel!("Time (days)", font=12)
  ylabel!("Abundance (nº individuals)", font=12)
  title!("H=0.1")
  savefig("SLC_NAt_H100_cij.png")
  
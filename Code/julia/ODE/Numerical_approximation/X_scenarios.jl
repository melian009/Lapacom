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
  t_0, k, r, K, H , d, Smax, gamma, cij = p
   
   
   
   #reproductive cycle
   #phi(t) = 2*pi*(t - t_0)/(365)
   #periodX(t) = 1/2*(1+tanh(2*sin(phi(t)) - k))
   
   function periodX(t)
    if (t % 365) / 365 >= 0.42
      return 1.0 # Reproductive Cycle
    else
      return 0.0 # Exploitation Cycle
    end
  end

   # reproductive capacity
   avg_size = du[2]
   Smat = 1.34 * (avg_size) - 28.06
   R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax - Smat)), 0.0), 1.0)

  du[1] = dNa = r * R_ * Na * ((K - Na)/ K) - d * Na - (1 - periodX(t)) * H * Na - cij * Na
  du[2] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - periodX(t))))
end





avg_oocytes = mean([92098, 804183]) # This is the actual mean. mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r_ = reggs*0.998611*0.971057*0.4820525*0.00629
# natural death rates per life stage.
d = 0.000322/365.14  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = 0.00014898749263737575
t0_ = 365.14*0.42
k_= 0.42
Naj_ = 2500
K_ = 64000.00    # Carrying capacity
Smax_ = 53.0             # Maximum size for adults3
#         t_0, k,  r,  K,  H ,  d, Smax,  gamma


n=10    #Number of years in the simulation
t_span = length(zeros(Float64,size(1:365.14*n)))
h_span = length(zeros(Float64, size(0:0.1:1)))
span = ones(Float64,size(1:365.14*n))
Kspan = ones(Float64,size(1:365.14*n))*K_ 
Sm_span = ones(Float64,size(1:365.14*n))*Smax_/2   # Linea de Smax/2 


Nai_cij = zeros(t_span,h_span)
Sai_cij = zeros(t_span,h_span)
periodX = zeros(t_span)


H_r = range(0, 1, length=h_span)
H_span = ones(Float64,h_span)
cij_span = ones(Float64,h_span)
N_H_c_span = zeros(11,11)
S_H_c_span = zeros(11,11)

for i in 1:length(H_span)
  H_span[i] = H_r[i]
  cij_span[i] = H_r[i]
end

Hs=length(zeros(Float64,size(1:length(H_span))))
cijs = length(zeros(Float64,size(1:length(cij_span))))


H1_ = H_span[1]
cij_= cij_span[1]

p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate,cij_,Naj_]  
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
    resultados_simulaciones = zeros(length(time2), length(H_span), length(cij_span)) 
    longitud_simulacion = length(time2)










for i in 1:length(H_span)
  H1_ = H_span[i]
  for j in 1:length(cij_span)

   cij_= cij_span[j]

    p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate,cij_,Naj_]  
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

    if i == 1 && j == 1
      longitud_simulacion = length(time2)
    # Guardar los resultados de la simulación en la matriz
      resultados_simulaciones[:,i,j] = vars[:,1]
    
    else 
      lon_lar = length(resultados_simulaciones[:,1,1])
      lon_cor = length(vars[:,1])
      elementos_faltantes = lon_lar - lon_cor
      if elementos_faltantes > 0
        vector_corto = vcat(vars[:,1],zeros(elementos_faltantes))
      else
        vector_corto = vars[:,1][1:lon_lar] # Trunca el vector corto si es más largo que el vector largo
      end
    resultados_simulaciones[:,i,j] = vector_corto
    end
  end
end


Na1c0=[resultados_simulaciones[:,1,1]  resultados_simulaciones[:,5,1]  resultados_simulaciones[:,11,1]]
Na5c0=[resultados_simulaciones[:,1,5],resultados_simulaciones[:,5,5], resultados_simulaciones[:,11,5]]
Na10c0=[resultados_simulaciones[:,1,10],resultados_simulaciones[:,5,10], resultados_simulaciones[:,11,10]]
Na11c0=[resultados_simulaciones[:,1,11],resultados_simulaciones[:,5,11], resultados_simulaciones[:,11,11]]

resultados_simulaciones[:,1,1]
length(resultados_simulaciones[:,1,5])


#X=0 en el dia 625 (day, H, cij)
X0 = resultados_simulaciones[625,:,1]
X1 = resultados_simulaciones[625,:,6]
X2 = resultados_simulaciones[625,:,10]
X3 = resultados_simulaciones[625,:,11]
#X=1 en el dia 900
X4 = resultados_simulaciones[900,:,1]
X5 = resultados_simulaciones[900,:,6]
X6 = resultados_simulaciones[900,:,10]
X7 = resultados_simulaciones[900,:,11]


plot(H_span,X0,label="cij=0.0", color=:blue)
plot!(H_span,X1,label="cij=0.5", color=:green)
plot!(H_span,X2,label="cij=0.9", color=:red)
plot!(H_span,X3,label="cij=1.0", color=:black,
      background=nothing)

xlabel!("Exploitation rate (H)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("X=0")
savefig("SLC_NA_H_cij_X0.png")


plot(H_span,X4,label="cij=0.0", color=:blue, style=:solid)
plot!(H_span,X5,label="cij=0.5", color=:green,style=:solid)
plot!(H_span,X6,label="cij=0.9", color=:red, style=:solid)
plot!(H_span,X7,label="cij=1.0", color=:black, style=:solid,
    background =nothing)
xlabel!("Exploitation rate (H)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("X=1")
savefig("SLC_NA_H_cij_X1.png")


plot(resultados_simulaciones[:,1,1],label="H=0.0, cij=0.0",  color=:blue,style=:solid)
plot!(resultados_simulaciones[:,1,6],label="H=0.5, cij=0.0",  color=:blue,style=:dash)
plot!(resultados_simulaciones[:,1,10],label="H=0.9, cij=0.0", color=:blue,style=:dashdot)

plot!(resultados_simulaciones[:,6,1],label="H=0.0, cij=0.5", color=:green,style=:solid)
plot!(resultados_simulaciones[:,6,6],label="H=0.5, cij=0.5", color=:green,style=:dash)
plot!(resultados_simulaciones[:,6,10],label="H=0.9, cij=0.5", color=:green,style=:dashdot)

plot!(resultados_simulaciones[:,11,1],label="H=0.0, cij=1.0", color=:red,style=:solid)
plot!(resultados_simulaciones[:,11,6],label="H=0.5, cij=1.0", color=:red,style=:dash)
plot!(resultados_simulaciones[:,11,10],label="H=0.9, cij=1.0", color=:red,style=:dashdot)

plot!(resultados_simulaciones[:,1,11],label="H=1, cij=0.0", color=:black,style=:solid)
plot!(resultados_simulaciones[:,6,11],label="H=1, cij=0.5", color=:black,style=:dash)
plot!(resultados_simulaciones[:,11,11],label="H=1,cij=1.0", color=:black,style=:dashdot,
       background=nothing,
       legend=:outerright)

xlims!(620,1000)
ylims!(0,7*10^4)
xlabel!("Time (days)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("X=[0,1]")
savefig("SLC_NA_H_cij_X.png")

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
    if (t % t_0 / t_0) >= k
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



avg_oocytes = [77404, 385613]
 # This is the actual mean.
reggs = avg_oocytes / (365 * 0.42) #aplication of the reproduction period stablish by law. (The time banned for extraction or exploitation for the species)
r_ = reggs*0.998611*0.971057*0.4820525*0.00629 # conversion rate of adults to eggs.
# natural death rates per life stage.
d = ([0.55,0.59])/365.14  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = [0.32,0.36] #0.00014898749263737575
t0_ = 365.14*0.42
k_= 0.42
Naj_ = 2500
K_ = 64000.00    # Carrying capacity
Smax_ = 53.0             # Maximum size for adults3
#         t_0, k,  r,  K,  H ,  d, Smax,  gamma


n=10    #Number of years in the simulation
t_span = length(zeros(Float64,size(1:365.14*n)))
h_span = length(zeros(Float64, size(0:0.01:1)))
span = ones(Float64,size(1:365.14*n))
#Kspan = ones(Float64,size(1:365.14*n))*K_ 
#Sm_span = ones(Float64,size(1:365.14*n))*Smax_/2   # Linea de Smax/2 


Nai_cij = zeros(t_span,h_span)
Sai_cij = zeros(t_span,h_span)
periodX = zeros(t_span)


H_r = range(0, 1, length=h_span)
H_span = ones(Float64,h_span)
cij_span = ones(Float64,h_span)
#N_H_c_span = zeros(11,11)
#S_H_c_span = zeros(11,11)

for i in 1:length(H_span)
  H_span[i] = H_r[i]
  cij_span[i] = H_r[i]
end

Hs=length(zeros(Float64,size(1:length(H_span))))
cijs = length(zeros(Float64,size(1:length(cij_span))))

#condiciones de estabilidad
H1_ = H_span[1]
cij_= cij_span[1]

n1_p_span = [t0_, k_, r_[1], K_, H1_, d[1], Smax_, size_growth_rate[1],cij_,Naj_]  
n=10  #Number of years in the simulation
tspan = (0,365.14*n)
U0_ = [10^4, 49.25]
prob_ = ODEProblem(SLC!, U0_, tspan, n1_p_span)
solve_= solve(prob_, Tsit5())

t_l = length(zeros(Float64,size(1:length(solve_.u))))
n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))
vars = zeros(t_l,n_v)
time_inicial = zeros(t_l)

for j in 1:length(solve_.u)
  for i in 1:2
    vars[j,i] = solve_.u[j][i]
    time_inicial[j] = solve_.t[j]
  end 
end

#Generacion de metrices cúbicas de almacenamiento 
resultados_simulaciones = zeros(length(time_inicial), length(H_span), length(cij_span)) 
tiempos_totales = zeros(length(time_inicial), length(H_span), length(cij_span))
tiempos_maximos = zeros(length(H_span), length(cij_span))


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
    time2 = zeros(t_l)
 
      for a in 1:length(solve_.u)
        for b in 1:2
          vars[a,b] = solve_.u[a][b]
          time2[a] = solve_.t[a]
        end 
      end

    if i == 1 && j == 1
      longitud_simulacion = length(time2)
    # Guardar los resultados de la simulación en la matriz cúbica
      resultados_simulaciones[:,i,j] = vars[:,1] # u = vars[:,1] = Abundancias // vars[:,2] = Tallas
      tiempos_totales[:,i,j] = time2
    else 
      lon_0 = length(resultados_simulaciones[:,1,1])
      lon_c = length(vars[:,1])
      elementos_faltantes = lon_0 - lon_c
      if elementos_faltantes > 0 # Concadena el vector corto con zeros para tener la misma longitud del vector de las condiciones iniciales, el mas largo.
        vector_corto = vcat(vars[:,1],zeros(elementos_faltantes))
        vector_corto_t = vcat(time2,ones(elementos_faltantes)*maximum(time2))
      else
        vector_corto = vars[:,1][1:lon_0] # Selecciona los valores del vector mas largo hasta la longitud de la simulación de las condiciones estandar (H=0,c=0)
        vector_corto_t = time2[1:lon_0]
      end
    resultados_simulaciones[:,i,j] = vector_corto
    tiempos_totales[:,i,j] = vector_corto_t
    tiempos_maximos[i,j] = maximum(vector_corto_t) #ultima iteración (nº máximo de días) por simulación
    end
  end
end

resultados_simulaciones
tiempos_totales
tiempos_maximos

h_0_c0_ = hcat(tiempos_totales[:,1,1],resultados_simulaciones[:,1,1]) #H =0, c=0
h_1_c0_ = hcat(tiempos_totales[:,101,1],resultados_simulaciones[:,101,1]) #H=0.99, c=0
h_0_c1_ = hcat(tiempos_totales[:,1,101],resultados_simulaciones[:,1,101]) #H=0, c=1
h_1_c1_ = hcat(tiempos_totales[:,101,101],resultados_simulaciones[:,101,101]) #H=0.99, C=1 

plot(h_0_c0_[:,1],h_0_c0_[:,2],label="H=0, c=0",  color=:blue, style=:solid)
plot!(h_1_c0_[:,1],h_1_c0_[:,2],label="H=1, c=0",  color=:red, style=:solid)
plot!(h_0_c1_[:,1],h_0_c1_[:,2],label="H=0, c=1",  color=:green, style=:solid)
plot!(h_1_c1_[:,1],h_1_c1_[:,2],label="H=1, c=1",  color=:brown, style=:solid)

xlims!(000,2000)
xlabel!("Time (days)", font=12)
ylabel!("Abundance (nº individuals)", font=12)

N1_mean(h_0_c0_[2])
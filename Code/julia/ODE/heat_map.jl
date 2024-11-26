using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics
using Random
using Distributions

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
  Na1, Na2, Sa1, Sa2 = u
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
   avg_size_1 = du[3]
   Smat_1 = 1.34 * (avg_size_1) - 28.06
   R_1 = min(max(0.5 * (1.0 + (avg_size_1 - Smat_1) / (Smax - Smat_1)), 0.0), 1.0)
   
   avg_size_2 = du[4]
   Smat_2 = 1.34 * (avg_size_1) - 28.06
   R_2 = min(max(0.5 * (1.0 + (avg_size_2 - Smat_2) / (Smax - Smat_2)), 0.0), 1.0)
  
  du[1] = dNa1 = r[1] * R_1 * Na1 * ((K - Na1)/ K) - d[1] * Na1 - (1 - periodX(t)) * H * Na1 - cij * Na2 #N1 
  du[2] = dNa2 = r[2] * R_2 * Na2 * ((K - Na2)/ K) - d[2] * Na2 - (1 - periodX(t)) * H * Na2 - cij * Na1 #N2
  du[3] = dSa1 = gamma[1] * Sa1 * (1 - Sa1 / (Smax - Smax * H * (1 - periodX(t))))
  du[4] = dSa2 = gamma[2] * Sa2 * (1 - Sa2 / (Smax - Smax * H * (1 - periodX(t))))

end



avg_oocytes = [77404, 385613] # This is the actual mean.
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
H1_ = H_span[101]
cij_= cij_span[101]

n1_p_span = [t0_, k_, r_[1], K_, H1_, d[1], Smax_, size_growth_rate[1],cij_]
n2_p_span = [t0_, k_, r_[2], K_, H1_, d[2], Smax_, size_growth_rate[2],cij_]  

p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate,cij_]



n=10  #Number of years in the simulation
tspan = (0,365.14*n)
U0_ = [10^4,10^4, mean([33.4,37.4]), mean([34.6,37.5])]
prob_ = ODEProblem(SLC!, U0_, tspan, p_span)
solve_= solve(prob_, Tsit5())


#Push variables in sigle vectors
N1_values = Float64[]
for i in 1:size(solve_.u , 1)
  push!(N1_values, solve_.u[i][1])  
end

N2_values = Float64[]
for i in 1:size(solve_.u , 1)
  push!(N2_values, solve_.u[i][2]) 
end
t_values = Float64[]
for i in 1:size(solve_.t , 1)
  push!(t_values, solve_.t[i][1])  
end


Random.seed(size_growth_rate[1]) # Setting the seed
d = Normal(μ=0.32, σ=0.05)  #Media, Varianza
n=rand(d,1000)             #Replicas aleatorioas, en este caso 1000 replicas.

N1_H1_C1  = [min(N1_values),max() # min, max
N2_H1_C1 # = mean(N2_values)

N1_H1_C0 = 0 # = mean(N1_values) # min, max
N2_H1_C0 = 0 # = mean(N2_values)


N1_H0_C0 # = mean(N1_values) # min, max
N2_H0_C0 # = mean(N2_values)
Mat_freq_rep

N1_H099_C0 # = mean(N1_values)
N2_H099_C0 # = mean(N2_values)

N1_H0_C099 #  = mean(N1_values)
N2_H0_C099 # = mean(N2_values)

N1_H099_C099 = mean(N1_values)
N2_H099_C099 = mean(N2_values)


heatmap()

#=
t_l = length(solve_.t)
n_v = length(solve_.u[1]) # Number of variables from the output of the ODE Problem Solve

#=
plot(t_values,N1_values)
plot!(t_values,N2_values)
ylims!(6.39*10^4,6.4*10^4)
t_l = length(zeros(Float64,size(1:length(solve_.u))))
n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))
=#

vars = zeros(t_l,n_v-2)  # Takes only abundances from both species
time_inicial = zeros(t_l)

for j in 1:length(N1_values)
  for i in 1:2
    vars[j,i] = solve_.u[j][i]
    time_inicial[j] = solve_.t[j]
  end 
end

#Generacion de metrices cúbicas de almacenamiento 
resultados_simulacionesN1 = zeros(length(time_inicial), length(H_span), length(cij_span)) #Matrix de salida para N1
resultados_simulacionesN2 = zeros(length(time_inicial), length(H_span), length(cij_span)) #Matrix de salida para N2

tiempos_totales = zeros(length(time_inicial), length(H_span), length(cij_span))
tiempos_maximos = zeros(length(H_span), length(cij_span),length(vars))



#=for i in 1:length(H_span)
  H1_ = H_span[1]
# for j in 1:length(cij_span)

   cij_= cij_span[2]
    
   p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate,cij_]  
   n=10  #Number of years in the simulation
   tspan = (0,365.14*n)
   U0_ = [10^4,10^4, mean([33.4,37.4]), mean([34.6,37.5])]
   prob_ = ODEProblem(SLC!, U0_, tspan, p_span)
   solve_= solve(prob_, Tsit5())

   t_l = length(zeros(Float64,size(1:length(solve_.u))))
    
   n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))


   vars = zeros(t_l,n_v-2)
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
      resultados_simulacionesN1[:,i,j] = vars[:,1] # N1: u = vars[:,1] = Abundancias // vars[:,3] = Tallas
      resultados_simulacionesN2[:,i,j] = vars[:,2] # N2: u = vars[:,2] = Abundancias // vars[:,4] = Tallas

      tiempos_totales[:,i,j] = time2
    else 
      lon_0_1 = length(resultados_simulacionesN1[:,1,1])
      lon_c_1 = length(vars[:,1])
      elementos_faltantes_1 = lon_0_1 - lon_c_1


      lon_0_2 = length(resultados_simulacionesN2[:,1,1])
      lon_c_2 = length(vars[:,2])
      elementos_faltantes_2 = lon_0_2 - lon_c_2

      if elementos_faltantes_1 > 0 || elementos_faltantes_2 > 0 # Concadena el vector corto con zeros para tener la misma longitud del vector de las condiciones iniciales, el mas largo.
        vector_corto_1 = vcat(vars[:,1],zeros(elementos_faltantes_1))
        vector_corto_t_1 = vcat(time2,ones(elementos_faltantes_1)*maximum(time2))
  
        vector_corto_2 = vcat(vars[:,2],zeros(elementos_faltantes_2))
        vector_corto_t_2 = vcat(time2,ones(elementos_faltantes_2)*maximum(time2))
      else
        vector_corto_1 = vars[:,1][1:lon_0_1] # Selecciona los valores del vector mas largo hasta la longitud de la simulación de las condiciones estandar (H=0,c=0)
        vector_corto_t_1 = time2[1:lon_0_1]

        vector_corto_2 = vars[:,2][1:lon_0_2] # Selecciona los valores del vector mas largo hasta la longitud de la simulación de las condiciones estandar (H=0,c=0)
        vector_corto_t_2 = time2[1:lon_0_2]
      end

     resultados_simulacionesN1[:,i,j] = vector_corto_1
     resultados_simulacionesN2[:,i,j] = vector_corto_2

    tiempos_totales[:,i,j] = vector_corto_t_1
    
    tiempos_maximos[i,j,1] = maximum(vector_corto_t_1) #ultima iteración (nº máximo de días) por simulación
    
    tiempos_maximos[i,j,2] = maximum(vector_corto_t_2) #ultima iteración (nº máximo de días) por simulación
    end
  end
end
=#
resultados_simulacionesN1[:,2,1]
resultados_simulacionesN2

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



using Random, Distributions
Random.seed(size_growth_rate[1]) # Setting the seed
d = Normal(μ=0.32, σ=0.05)  #Media, Varianza
n=rand(d,1000)             #Replicas aleatorioas, en este caso 1000 replicas.

gr()
data = rand(21,100)
heatmap(1:size(data,1),
    1:size(data,2), data,
    c=cgrad([:blue, :white,:red, :yellow]),
    xlabel="x values", ylabel="y values",
    title="My title")

= #
using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics
using Random
using Distributions
using StatsPlots



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
   # N_1
   avg_size_1 = du[3]
   Smat_1 = 1.34 * (avg_size_1) - 28.06
   R_1 = min(max(0.5 * (1.0 + (avg_size_1 - Smat_1) / (Smax - Smat_1)), 0.0), 1.0)
   # N_2
   avg_size_2 = du[4]
   Smat_2 = 1.34 * (avg_size_1) - 28.06
   R_2 = min(max(0.5 * (1.0 + (avg_size_2 - Smat_2) / (Smax - Smat_2)), 0.0), 1.0)
  
  du[1] = dNa1 = r[1] * R_1 * Na1 * ((K - Na1)/ K) - d[1] * Na1 - (1 - periodX(t)) * H * Na1 - cij * Na2 * Na1 #N1 
  du[2] = dNa2 = r[2] * R_2 * Na2 * ((K - Na2)/ K) - d[2] * Na2 - (1 - periodX(t)) * H * Na2 - cij * Na1 * Na2 #N2
  du[3] = dSa1 = gamma[1] * Sa1 * (1 - Sa1 / (Smax - Smax * H * (1 - periodX(t))))
  du[4] = dSa2 = gamma[2] * Sa2 * (1 - Sa2 / (Smax - Smax * H * (1 - periodX(t))))

end




avg_oocytes = [77404, 385613] # This is the actual mean.
reggs = avg_oocytes / (365 * 0.42) #aplication of the reproduction period stablish by law. (The time banned for extraction or exploitation for the species)
r_ = reggs*0.998611*0.971057*0.4820525*0.00629 # conversion rate of adults to eggs.
# natural death rates per life stage.
d_ = ([0.55,0.59])/365.14  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = [0.32,0.36] #0.00014898749263737575
t0_ = 365.14*0.42
k_= 0.42

K_ = 64000.00    # Carrying capacity
Smax_ = 53.0             # Maximum size for adults
#         t_0, k,  r,  K,  H ,  d, Smax,  gamma


n_simulaciones = 10 # Número de simulaciones
t_span = (0.0, 365.14*10)  # Tiempo de simulación (por ejemplo, un año)
t_plt = 0.0:1.0:365.14*10  # Los tiempos en los que se evaluará la solución


h_span = length(zeros(Float64, size(0:0.1:1)))
H_r = range(0, 1, length=h_span)
H_span = ones(Float64,h_span)
cij_span = ones(Float64,h_span)

for i in 1:length(H_span)
  H_span[i] = H_r[i]
  cij_span[i] = H_r[i]
end




j=1
n=11


# for j in 1:length(H_span)
  H = H_span[j] #Exploitation
#  for k in 1:length(cij_span)
  cij = cij_span[n]  #Simetric competence component

    
    # Almacenar los conjuntos resultados de las simulaciones
     resultados_t_concatenados = Float64[]  # Para almacenar los valores de t
    
     resultados_Na1_concatenados = Float64[]  # Para almacenar los valores de Na1
     resultados_Na2_concatenados = Float64[]  # Para almacenar los valores de Na2
     resultados_Sa1_concatenados = Float64[]  # Para almacenar los valores de Sa1
     resultados_Sa2_concatenados = Float64[]  # Para almacenar los valores de Sa2

    # # Almacenar los conjuntos resultados de las simulaciones
     resultados_t = Float64[]  # Para almacenar los valores de t
    
    resultados_Na1 = Float64[]  # Para almacenar los valores de Na1
     resultados_Na2 = Float64[]  # Para almacenar los valores de Na2
     resultados_Sa1 = Float64[]  # Para almacenar los valores de Sa1
     resultados_Sa2 = Float64[]  # Para almacenar los valores de Sa2

    # Realizamos las simulaciones
     for i in 1:n_simulaciones
      # Generamos valores aleatorios para los parámetros (distribución normal)
      t_0 = t0_ + 0.0001 * randn()
      k = k_ + 0.01 * randn()
      r = [r_[1] + 0.01 * randn(), r_[2] + 0.01 * randn()] 
      K = K_ + 0.1 * randn()
      gamma = [size_growth_rate[1] + 0.01 * randn(),size_growth_rate[2] + 0.01 * randn()]
      d = [d_[1], d_[2]]
      Smax = Smax_ + 0.1 * randn()
    
      # Condiciones iniciales
      U0_ = [10^4,10^4, mean([33.4,37.4]), mean([34.6,37.5])]
    
      # Definir el problema diferencial
      prob = ODEProblem(SLC!, U0_, t_span, [t_0, k, r, K, H, d, Smax, gamma, cij])
    
      #  Resolver el problema
      sol = solve(prob, Tsit5())
    
      # Almacenar las soluciones de t, Na1, Na2, Sa1, Sa2
      for m in 1:size(sol.t , 1)
        push!(resultados_t, sol.t[m])
        push!(resultados_Na1, sol.u[m][1])
        push!(resultados_Na2, sol.u[m][2])
        push!(resultados_Sa1, sol.u[m][3])
        push!(resultados_Sa2, sol.u[m][4])
      end
    
      # Concatenar los resultados para obtener las soluciones completas de Na1, Na2, Sa1, Sa2

         resultados_t_concatenados = vcat(resultados_t...)
    
         resultados_Na1_concatenados = vcat(resultados_Na1...)
         resultados_Na2_concatenados = vcat(resultados_Na2...)
         resultados_Sa1_concatenados = vcat(resultados_Sa1...)
         resultados_Sa2_concatenados = vcat(resultados_Sa2...)
    end 

     
    plot(resultados_Na1_concatenados, resultados_Na2_concatenados, xlabel= "N1", ylabel = "N2", label=vcat("H=", H, "CIJ=",cij))
    
    plot!(resultados_Na1_concatenados, resultados_Na2_concatenados, label=vcat("H=", H, "CIJ=",cij))
    
#  end
#end

plot(resultados_t,log.(resultados_Na1), xlabel= "t", ylabel = "N1")

# Generar la distribución de frecuencias de Na1 y Na2
#Populations
density(log.(resultados_Na1_concatenados), label="Na1", ylabel="Frecuencia", title="Distribución de frecuencias de Na1")
density!(log.(resultados_Na2_concatenados), label="Na2", ylabel="Frecuencia", title="Distribución de frecuencias de Na2")
xlims!(6.0*10^4,7.0*10^4)

#Size
density(log.(resultados_Sa1_concatenados), label="Na1", ylabel="Frecuencia", title="Distribución de frecuencias de Na2")
density!(log.(resultados_Sa2_concatenados), label="Na2", ylabel="Frecuencia", title="Distribución de frecuencias de Na2")
xlims!(3.9,3.98)


# Obtener el mínimo y máximo de las soluciones de las ecuaciones
min_max_Na1 = [minimum(resultados_Na1_concatenados), maximum(resultados_Na1_concatenados)]
min_max_Na2 = [minimum(resultados_Na2_concatenados), maximum(resultados_Na2_concatenados)]
min_max_Sa1 = [minimum(resultados_Sa1_concatenados), maximum(resultados_Sa1_concatenados)]
min_max_Sa2 = [minimum(resultados_Sa2_concatenados), maximum(resultados_Sa2_concatenados)]

# Mostrar los resultados mínimos y máximos
println("Mínimo y máximo de Na1: ", min_max_Na1)
println("Mínimo y máximo de Na2: ", min_max_Na2)
println("Mínimo y máximo de Sa1: ", min_max_Sa1)
println("Mínimo y máximo de Sa2: ", min_max_Sa2)

density(resultados_Na1_concatenados, bins=300, label="Na1", ylabel="Frecuencia", title="Distribución de frecuencias de Na1")
density!(resultados_Na2_concatenados, bins=300, label="Na2", ylabel="Frecuencia", title="Distribución de frecuencias de Na2")
xlims!(min_max_Na1[2]-10^3,min_max_Na2[2]+100)


# Rango y la Resolución del grid para Na1 y Na2
  n_bins = 100  # Número de bins para el histograma (resolución)
  x_min, x_max = minimum(resultados_Na1_concatenados), maximum(resultados_Na1_concatenados)
  y_min, y_max = minimum(resultados_Na2_concatenados), maximum(resultados_Na2_concatenados)
  
  # Grid de bins para Na1 y Na2
  x_bins = range(x_min, stop=x_max, length=n_bins)
  y_bins = range(y_min, stop=y_max, length=n_bins)
  
  # Matriz para almacenar las frecuencias
  frequencies = zeros(Int, n_bins, n_bins)
  
  # Contabilización de las frecuencias de cada pares (Na1, Na2)
  for i in 1:length(resultados_Na1_concatenados)
      # Localización del bin correspondiente para Na1 y Na2
      x_bin = searchsortedlast(x_bins, resultados_Na1_concatenados[i])
      y_bin = searchsortedlast(y_bins, resultados_Na2_concatenados[i])
  
      # Incremento del contador en la posición correspondiente
      frequencies[x_bin, y_bin] += 1
  end
  
  # Matriz de frecuencias normalizada 
  frequencies_norm = frequencies / sum(frequencies)
  
  # Heatmap
  heatmap(x_bins, y_bins, frequencies_norm, xlabel="Na1", ylabel="Na2", title="Heatmap de Frecuencia (Na1 vs Na2)", color=:viridis, clims=(0.8, 1))
  savefig("Heatmap_h0_c0.png")
  # Heatmap LOG SCALE
  heatmap(x_bins, y_bins, log.(frequencies_norm.^(-1)), xlabel="Na1", ylabel="Na2", title="Heatmap de Frecuencia (Na1 vs Na2)", color=:viridis)
  savefig("Heatmap_log_h0_c0.png")
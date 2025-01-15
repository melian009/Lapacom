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
Cij_r = range(0, 1, length=h_span)

H_span = ones(Float64,h_span)
cij_span = ones(Float64,h_span)

for i in 1:length(H_span)
  H_span[i] = H_r[i]
  cij_span[i] = Cij_r[i]
end


#Fig 4a: with discrete leyend.

for j in 1:11
  for n in 1:11 
  H = H_span[j] #Exploitation

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
      sol = solve(prob, maxiters=500)
    
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

    if j == 1 && n == 1 
      plot(resultados_Na1_concatenados, resultados_Na2_concatenados, xlabel= "N1", ylabel = "N2", label=vcat("H=", H, "CIJ=",cij), legend=:outerright)
    else 
      plot!(resultados_Na1_concatenados, resultados_Na2_concatenados, label=vcat("H=", H, "CIJ=",cij))
    end
end
end
plot!(legend=:none)


#Heatmaps: by specific conditions

j = 1
n = 1 

H = H_span[j] #Exploitation
cij = cij_span[n]  #Simetric competence component

    
# Almacenar los conjuntos resultados de las simulaciones
resultados_t_concatenados = Float64[]  # Para almacenar los valores de t
    
resultados_Na1_concatenados = Float64[]  # Para almacenar los valores de Na1
resultados_Na2_concatenados = Float64[]  # Para almacenar los valores de Na2
resultados_Sa1_concatenados = Float64[]  # Para almacenar los valores de Sa1
resultados_Sa2_concatenados = Float64[]  # Para almacenar los valores de Sa2

# Almacenar los conjuntos resultados de las simulaciones
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
  sol = solve(prob, maxiters=500)
    
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

# Rango y la Resolución del grid para Na1 y Na2
n_bins = 100 # Número de bins para el histograma (resolución)
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
Bin_pos = argmax(frequencies_norm)

#= if H==1 && cij==1
 anotations_00 = hcat(cij,H,x_bins[Bin_pos[1]], y_bins[Bin_pos[2]],frequencies_norm[Bin_pos])
 MATRIX = (anotations_00)
 # else
  anotations_ = hcat(cij,H,x_bins[Bin_pos[1]], y_bins[Bin_pos[2]],frequencies_norm[Bin_pos])
  MATRIX = vcat(MATRIX,anotations_)
 # end =# 

# Heatmap
 heatmap(x_bins,
 y_bins, 
 frequencies_norm,
 xlabel="Na1", 
 ylabel="Na2", 
 title="(Na1 vs Na2)", 
 color=cgrad(:viridis, rev=true),
 clims=(minimum(frequencies_norm), maximum(frequencies_norm)))
 #xlims!(0,6.4*10^4)
 #ylims!(0,6.4*10^4)

 #annotate!(x_bins[Bin_pos[1]], y_bins[Bin_pos[2]], text("H=$H", 8))

# savefig("Heatmap_h0_4_c0.png")
  
# Heatmap LOG (SCALE^-1)
 heatmap(x_bins,
 y_bins, 
 log.(frequencies_norm), 
 xlabel="Na1",
 ylabel="Na2",
 title=vcat("H=", H, "CIJ=",cij," (Na1 vs Na2)"),
 color=:viridis,
 clims=(minimum(log.(frequencies_norm.^(-1))), maximum(log.(frequencies_norm.^(-1)))))
  
# savefig("Heatmap_log_h0_c0.png")






#Density diagrams for limit cycles
#Frecuency distribution of: 
#Populations Na1 y Na2
density((resultados_Na1_concatenados), label="Na1", ylabel="Frecuencia")
density!((resultados_Na2_concatenados), label="Na2", ylabel="Frecuencia")
xlims!(0,7.0*10^4)

#Size Sa1 y Sa2

density((resultados_Sa1_concatenados), label="Na1", ylabel="Frecuencia", title="Distribución de frecuencias de Na2")
density!((resultados_Sa2_concatenados), label="Na2", ylabel="Frecuencia", title="Distribución de frecuencias de Na2")
xlims!(3.9,3.98)

#Obtain min and max of the equation solutions
min_max_Na1 = [minimum(resultados_Na1_concatenados), maximum(resultados_Na1_concatenados)]
min_max_Na2 = [minimum(resultados_Na2_concatenados), maximum(resultados_Na2_concatenados)]
min_max_Sa1 = [minimum(resultados_Sa1_concatenados), maximum(resultados_Sa1_concatenados)]
min_max_Sa2 = [minimum(resultados_Sa2_concatenados), maximum(resultados_Sa2_concatenados)]

#Showe min and max outputs
println("Mínimo y máximo de Na1: ", min_max_Na1)
println("Mínimo y máximo de Na2: ", min_max_Na2)
println("Mínimo y máximo de Sa1: ", min_max_Sa1)
println("Mínimo y máximo de Sa2: ", min_max_Sa2)

#Frecuency distribution of Populations Na1 y Na2 with title
density(resultados_Na1_concatenados, bins=300, label="Na1", ylabel="Frecuencia", title="Distribución de frecuencias de Na1")
density!(resultados_Na2_concatenados, bins=300, label="Na2", ylabel="Frecuencia", title="Distribución de frecuencias de Na2")





# Heatmaps consolidados
# Variables para almacenar las frecuencias consolidadas
n_bins = 100
consolidated_frequencies_H = zeros(Float64, n_bins, n_bins)
consolidated_frequencies_cij = zeros(Float64, n_bins, n_bins)

surf_plots = []
heat_plots = []
j=1
n=1

#for j in 1:length(cij_span)
    for n in 1:length(H_span)
        H = H_span[n]  # Valor de explotación
        cij = cij_span[j]  # Componente de competencia simétrica
        
        # Reiniciar resultados para esta combinación de H y cij
        resultados_Na1_concatenados =Float64[]
        resultados_Na2_concatenados =Float64[]

        for i in 1:n_simulaciones
            # Generar valores aleatorios para parámetros
            t_0 = t0_ + 0.0001 * randn()
            k = k_ + 0.01 * randn()
            r = [r_[1] + 0.01 * randn(), r_[2] + 0.01 * randn()] 
            K = K_ + 0.1 * randn()
            gamma = [size_growth_rate[1] + 0.01 * randn(), size_growth_rate[2] + 0.01 * randn()]
            d = [d_[1], d_[2]]
            Smax = Smax_ + 0.1 * randn()

            # Condiciones iniciales
            U0_ = [10^4, 10^4, mean([33.4, 37.4]), mean([34.6, 37.5])]

            # Resolver el problema diferencial
            prob = ODEProblem(SLC!, U0_, t_span, [t_0, k, r, K, H, d, Smax, gamma, cij])
            sol = solve(prob, maxiters=500)

            # Almacenar resultados de Na1 y Na2
            for m in 1:length(sol.t)
                push!(resultados_Na1_concatenados, sol.u[m][1])
                push!(resultados_Na2_concatenados, sol.u[m][2])
            end
        end
        # Rango y la Resolución del grid para Na1 y Na2
        x_min, x_max = minimum(resultados_Na1_concatenados), maximum(resultados_Na1_concatenados)
        y_min, y_max = minimum(resultados_Na2_concatenados), maximum(resultados_Na2_concatenados)
        # Grid de bins para Na1 y Na2
        x_bins = range(x_min, stop=x_max, length=n_bins)
        y_bins = range(y_min, stop=y_max, length=n_bins)
        # Contabilizar frecuencias para esta combinación de H y cij
        local_frequencies = zeros(Int, n_bins, n_bins)
        # Grid de bins para Na1 y Na2
  
        for i in 1:length(resultados_Na1_concatenados)
          # Localización del bin correspondiente para Na1 y Na2
          x_bin = searchsortedlast(x_bins, resultados_Na1_concatenados[i])
          y_bin = searchsortedlast(y_bins, resultados_Na2_concatenados[i])
      
          # Verificar que los índices estén dentro de los límites
          if x_bin > 0 && x_bin <= n_bins && y_bin > 0 && y_bin <= n_bins
              local_frequencies[x_bin, y_bin] += 1
          end
      end

        # Normalizar frecuencias locales
        local_frequencies_norm = local_frequencies / sum(local_frequencies)
        local_frequencies_norm

        # Acumular en las frecuencias consolidadas
        consolidated_frequencies_H .+= local_frequencies_norm
        consolidated_frequencies_cij .+= local_frequencies_norm
      # Crear gráfico de superficie
    if n == 1
      surf_0 = surface(
          x_bins, y_bins, consolidated_frequencies_H,
          xlabel = "Na1", ylabel = "Na2", 
          title = "P(Na(H))", 
          color = cgrad(:viridis, rev=true),
          clims = (minimum(consolidated_frequencies_H), maximum(consolidated_frequencies_H))
      )
    else
      surface!(
          x_bins, y_bins, consolidated_frequencies_H,
          xlabel = "Na1", ylabel = "Na2", 
          title = "P(Na(H))", 
          color = cgrad(:viridis, rev=true),
          clims = (minimum(consolidated_frequencies_H), maximum(consolidated_frequencies_H))
      )
    end
    if j == 1
      heat_0 = heatmap(
        x_bins, y_bins, consolidated_frequencies_cij,
        xlabel = "Na1", ylabel = "Na2",
        title = "P(Na(Cij))", 
        color = cgrad(:viridis),
        clims = (minimum(consolidated_frequencies_cij), maximum(consolidated_frequencies_cij))
        )
    else
      heat_0 = heatmap!(
      x_bins, y_bins, consolidated_frequencies_cij,
      xlabel = "Na1", ylabel = "Na2",
      title = "P(Na(Cij))", 
      color = cgrad(:viridis),
      clims = (minimum(consolidated_frequencies_cij), maximum(consolidated_frequencies_cij))
      )
    end
  end

#end


# Graficar heatmaps consolidados
# P([dNA]
surf_0
heat_0


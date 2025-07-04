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

n_simulaciones = 100 # Número de simulaciones
t_span = (0.0, 365.14*10)  # Tiempo de simulación (por ejemplo, un año)
t_plt = 0.0:1.0:365.14*10  # Los tiempos en los que se evaluará la solución (10 years)


h_span = length(0:0.1:1)  # Número de elementos en el rango
H_r = range(0, 1, length=h_span)  # Rango de H
Cij_r = range(0, 1, length=h_span)  # Rango de Cij


# Convertir los rangos en vectores
H_span = collect(H_r)
cij_span = collect(Cij_r)
 
limt_cycle = plot()

#Fig 4a: with discrete leyend
for j in 1:11 # Cij
cij = cij_span[j]  #Simetric competence component

  for n in 1:11 # H 
  H = H_span[n] #Exploitation
    
    
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
      t_0 = t0_ + t0_/3 * randn()
      k = k_ + k_/3 * randn()
      r = [r_[1] + r_[1]/3 * randn(), r_[2] + r_[2]/3 * randn()] 
      K = K_ + K_/30 * randn()
      gamma = [size_growth_rate[1] + size_growth_rate[1]/3 * randn(),size_growth_rate[2] + size_growth_rate[2]/3 * randn()]
      d = [d_[1]+d_[1]/3*randn(), d_[2]+d_[2]/3*randn()]
      Smax = Smax_ + Smax_/3 * randn()
    
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
    # Graficar los resultados
    if j == 1 && n == 1
      limt_cycle = plot!(resultados_Na1_concatenados, resultados_Na2_concatenados, 
             xlabel="N1", ylabel="N2", 
             label="H=$H, Cij=$cij",  color = cgrad(:viridis,H),linewidth=5)
    else
      limt_cycle =  plot!(resultados_Na1_concatenados, resultados_Na2_concatenados, 
              label="H=$H, Cij=$cij", linewidth=5)
    end
    display(limt_cycle)
  end
end

plot!(legend=:outerright)

xlims!(0, 7.0*10^4)
ylims!(0, 7.0*10^4)
plot!(background_color=:transparent, grid=true)
savefig("Figure_4a.png")


# Almacenar los conjuntos resultados de las simulaciones
n_simulaciones = 10 # Número de simulaciones
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
i=1
j=1
cij = cij_span[j]  #Simetric competence component
H = H_span[i] #Exploitation
for i in 1:n_simulaciones
  # Generamos valores aleatorios para los parámetros (distribución normal)
  t_0 = t0_ + t0_/3 * randn()
  k = k_ + k_/3 * randn()
  r = [r_[1] + r_[1]/3 * randn(), r_[2] + r_[2]/3 * randn()] 
  K = K_ + K_/30 * randn()
  gamma = [size_growth_rate[1] + size_growth_rate[1]/3 * randn(),size_growth_rate[2] + size_growth_rate[2]/3 * randn()]
  d = [d_[1]+d_[1]/3*randn(), d_[2]+d_[2]/3*randn()]
  Smax = Smax_ + Smax_/3 * randn()
    
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
x_min_2, x_max_2 = minimum(resultados_Sa1_concatenados), maximum(resultados_Sa1_concatenados)
y_min_2, y_max_2 = minimum(resultados_Sa2_concatenados), maximum(resultados_Sa2_concatenados)
  
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
 title=vcat("H= $H, Cij= $cij","(Na1 vs Na2)"),
 color=:viridis,
 clims=(minimum(log.(frequencies_norm.^(-1))), maximum(log.(frequencies_norm.^(-1)))))
 savefig("Heatmap_log_h0_c0.png")






#Density diagrams for limit cycles
#Frecuency distribution of: 
#Populations Na1 y Na2
density((resultados_Na1_concatenados), label="Na1", xlabel="Abundance", ylabel="Frecuency")
density!((resultados_Na2_concatenados), label="Na2")
xlims!(10^4,7.0*10^4)

#Size Sa1 y Sa2

density((resultados_Sa1_concatenados), label="Sa Patella aspera", xlabel="Size", ylabel="Frecuency", title="Size frequency")
density!((resultados_Sa2_concatenados), label="Sa Patella ordinaria")
#xlims!(3.9,3.98)

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

#Frecuency distribution of Size Sa1 y sa2 
density(resultados_Sa1_concatenados.^(-1), bins=300, label="Sa1", ylabel="Frecuencia", title="Distribución de frecuencias de Na1")
density!(resultados_Sa2_concatenados.^(-1), bins=300, label="Sa2", ylabel="Frecuencia", title="Distribución de frecuencias de Na2")





# Heatmaps consolidados
# Variables para almacenar las frecuencias consolidadas
n_simulaciones = 100
n_bins = 100
consolidated_frequencies_ = zeros(Float64, n_bins, n_bins)
heat_0=heatmap()
j=1
#n=1

for j in 1
  cij = cij_span[j]  # Componente de competencia simétrica
  for n in 1:11
      H = H_span[n]  # Valor de explotación      
      # Reiniciar resultados para esta combinación de H y cij
      resultados_Na1_concatenados = Float64[]
      resultados_Na2_concatenados = Float64[]

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
          sol = solve(prob, maxiters=1000)

          # Almacenar resultados de Na1 y Na2
          for m in 1:length(sol.t)
              push!(resultados_Na1_concatenados, sol.u[m][1])
              push!(resultados_Na2_concatenados, sol.u[m][2])
          end
      end

      # Rango y la Resolución del grid para Na1 y Na2
      x_min, x_max = 0, 6.7*10^4
      y_min, y_max = 0, 6.7*10^4
      # Grid de bins para Na1 y Na2
      x_bins = range(x_min, stop=x_max, length=n_bins)
      y_bins = range(y_min, stop=y_max, length=n_bins)
      # Contabilizar frecuencias para esta combinación de H y cij
      local_frequencies = zeros(Int, n_bins, n_bins)

      for i in 1:length(resultados_Na1_concatenados)
          # Localización del bin correspondiente para Na1 y Na2
          x_bin = searchsortedlast(x_bins, resultados_Na1_concatenados[i])
          y_bin = searchsortedlast(y_bins, resultados_Na2_concatenados[i])
      
          # Verificar que los índices estén dentro de los límites
          if x_bin > 0 && x_bin <= n_bins && y_bin > 0 && y_bin <= n_bins
              local_frequencies[x_bin, y_bin] += 1
          end
      end

      # Normalizar frecuencias 
      frequencies_norm = local_frequencies / sum(local_frequencies)
      # Acumular en las frecuencias consolidadas
      consolidated_frequencies_.+= frequencies_norm
      # Crear heatmap
      if j == 1 && n == 1
        heat_0 = heatmap(x_bins,
         y_bins, 
         consolidated_frequencies_,
         xlabel="Na1", 
         ylabel="Na2", 
         title="(Na1 vs Na2)", 
         color=cgrad(:thermal, rev=true),
         clims=(minimum(consolidated_frequencies_), maximum(consolidated_frequencies_)))
      else
        heat_0 = heatmap!(x_bins,
         y_bins, 
         consolidated_frequencies_,
         xlabel="Na1", 
         ylabel="Na2", 
         title="(Na1 vs Na2)", 
         color=cgrad(:thermal, rev=true),
         clims=(minimum(consolidated_frequencies_), maximum(consolidated_frequencies_)))
      end
      display(heat_0)
  end
end

# Mostrar gráfico final
heat_0
surface!(background_color=:transparent, grid=true)
savefig("Figure_4b_surface.png")





N_span = 11  # Número de elementos en el rango
H_span = range(0, 1, length=N_span) |> collect
cij_span = range(0, 1, length=N_span) |> collect

# Matrices para almacenar la abundancia promediada de cada especie
abundance_matrix_N1 = zeros(length(H_span), length(cij_span))
abundance_matrix_N2 = zeros(length(H_span), length(cij_span))

for j in 1:N_span
    cij = cij_span[j]  # Componente de competencia simétrica

    for n in 1:N_span
        H = H_span[n]  # Valor de explotación
        
        # Variables para acumular abundancias
        total_N1 = 0.0
        total_N2 = 0.0
        total_count = 0

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
            sol = solve(prob, maxiters=1000)

            # Calcular la abundancia final (o promedio en el tiempo)
            Na1_final = mean([sol.u[m][1] for m in 1:length(sol.t)])
            Na2_final = mean([sol.u[m][2] for m in 1:length(sol.t)])

            # Sumar a las abundancias totales
            total_N1 += Na1_final
            total_N2 += Na2_final
            total_count += 1
        end

        # Promediar abundancia y almacenar en la matriz correspondiente
        abundance_matrix_N1[n, j] = total_N1 / total_count
        abundance_matrix_N2[n, j] = total_N2 / total_count
    end
end

# Escalar valores para mejorar la visibilidad en el heatmap
abundance_matrix_N1 .*= 1e+4  # Ajusta este factor según el orden de magnitud de los datos
abundance_matrix_N2 .*= 1e+4

# Límites de color (para mejorar la visualización)
clims_N1 = (minimum(abundance_matrix_N1[(2:11),(2:11)]), maximum(abundance_matrix_N1[2:11]))
clims_N2 = (minimum(abundance_matrix_N2[2:11]), maximum(abundance_matrix_N2[2:11]))

# Crear heatmaps con escala ajustada
p1 = heatmap(H_span, cij_span, abundance_matrix_N1,
             xlabel="H", ylabel="Cij", 
             title="Abundancia Promediada de N1",
             color=cgrad(:thermal, rev=true),
             clims=clims_N1)

p2 = heatmap(H_span, cij_span, abundance_matrix_N2,
             xlabel="H", ylabel="Cij", 
             title="Abundancia Promediada de N2",
             color=cgrad(:thermal, rev=true))

display(p1)
display(p2)

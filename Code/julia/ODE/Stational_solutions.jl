using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics
using Random
using Distributions
using StatsPlots


```
Parameters and variables:
 - Ne = eggs abundance
 - Nt = trocophore abuncance
 - Nv = veliger abuncance
 - Nj = juvenile abuncance
 - Na = adults abuncance

 - Sa = adults size (Average sizes Before MPS+FULL = [46.44,44.45])
 - r = population growth rate [9.17,5.03]
 - R = reproductive capacity
 - K = carrying capacity (k = 1e^4)
 - X = Reproductive period [1,0] 
 - (1-X) = Exploitation period [1,0]
 - H = Exploitation rate (H = [0.639,0.57])
 - gEA = instant conversion rate of life stages (EA = Eggs to Adults) (gEA = 0.006)
 - de = natural mortality rate or death rate for eggs (de = [de_po,de_pa]) # Note: this value needs to be defined.
 - da = natural mortality rate or death rate for adults (da = [0.55,0.59]) # Note: empirical estimated values
 - Smax = maximum adult size estimated (56.0mm) # Note: Empirical value from the sample analized
 - gamma = adult growth rate (gamma=[0.32,0.36] year^{-1})
```

#Simple Life Cycle ([Na1,Sa1], Patella aspera, [Na2,Sa2], Patella ordinaria)

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
     # N1 Patella aspera
     avg_size_1 = du[3]
     Smat_1 = 1.34 * (avg_size_1) - 28.06
     R_1 = min(max(0.5 * (1.0 + (avg_size_1 - Smat_1) / (Smax - Smat_1)), 0.0), 1.0)
     # N2 Patella Ordinaria
     avg_size_2 = du[4]
     Smat_2 = 1.34 * (avg_size_1) - 28.06
     R_2 = min(max(0.5 * (1.0 + (avg_size_2 - Smat_2) / (Smax - Smat_2)), 0.0), 1.0)
    
    du[1] = dNa1 = r[1] * R_1 * Na1 * ((K - Na1)/ K) - d[1] * Na1 - (1 - periodX(t)) * H * Na1 - cij * Na2 * Na1 #N1 
    du[2] = dNa2 = r[2] * R_2 * Na2 * ((K - Na2)/ K) - d[2] * Na2 - (1 - periodX(t)) * H * Na2 - cij * Na1 * Na2 #N2
    du[3] = dSa1 = gamma[1] * Sa1 * (1 - Sa1 / (Smax - Smax * H * (1 - periodX(t))))
    du[4] = dSa2 = gamma[2] * Sa2 * (1 - Sa2 / (Smax - Smax * H * (1 - periodX(t))))
  
end

```
Empirical estimated mortalities (Z,d,F) and exploitation rates (E):

 {Patella ordinaria} (Henriques, 2012)

 - Natural Mortality rate (d) was 0.55 year^{-1}
 - Fishing mortality rate (F) was 1.24 year^{−1} 
 - Total mortality rate (Z=d+F) was 1.79 year^{−1}; 
 - Exploitation rate (E=F/Z) was 0.693. 

 {Patella aspera} (Sousa, 2017)

 - Natural mortality (d) was 0.59 year^{-1};
 - Fishing mortality (F) was 0.79 year^{-1};
 - Total mortality (Z=d+F) was 1.38 year^{-1};
 - Exploitation rate (E=F/Z) was 0.57. 

For numerical simulation use these rates (d and E).

Average sizes before and after marine protected area implementations

 Before (1996-2006): FULL ACCESS.
 Patella apera = 43.53mm
 Patella ordinaria = 46.26mm

 After  (2007-2017): MPA+FULL
 Patella aspera = 44.45mm
 Patella ordinaria = 46.44mm

 Only MPS (2007-2017)
 Patella aspera = 50.61mm
 Patella ordinaria = 49.25mm

 Only Full Access (2007-2017) 

 Patella aspera = 43.41mm
 Patella ordinaria = 45.72mm

```
  
#   Parameters for SLC

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
t_plt = 0.0:1.0:365.14*10  # Los tiempos en los que se evaluará la solución


h_span = length(zeros(Float64, size(0:0.1:1)))
H_r = range(0, 1, length=h_span)
Cij_r = range(0, 1, length=h_span)

H_span = ones(Float64,h_span)
cij_span = ones(Float64,h_span)

for i in 1:length(h_span)
  H_span[i] = H_r[i]
  cij_span[i] = Cij_r[i]
end


#Fig 4a: with discrete leyend

 for j in 1:11 # Cij
    cij = cij_span[j]  #Simetric competence component
    
    for n in 1:11 # H 
        H = H_span[n] #Exploitation
        
        # Almacenar los conjuntos resultados de cada simulación
        
        resultados_t = Float64[]  # Para almacenar los valores de t
        
        resultados_Na1 = Float64[]  # Para almacenar los valores de Na1
        resultados_Na2 = Float64[]  # Para almacenar los valores de Na2

        resultados_Sa1 = Float64[]  # Para almacenar los valores de Sa1
        resultados_Sa2 = Float64[]  # Para almacenar los valores de Sa2

        # Almacenar los conjuntos resultados de todas las simulaciones del escenario
        resultados_t_concatenados = Float64[]  # Para almacenar los valores de t
        
        resultados_Na1_concatenados = Float64[]  # Para almacenar los valores de Na1
        resultados_Na2_concatenados = Float64[]  # Para almacenar los valores de Na2
        resultados_Sa1_concatenados = Float64[]  # Para almacenar los valores de Sa1
        resultados_Sa2_concatenados = Float64[]  # Para almacenar los valores de Sa2
    
       
    
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
        
        #Condiciones iniciales
        U0_ = [10^4,10^4, mean([33.4,37.4]), mean([34.6,37.5])]
        
        #Definir el problema diferencial
        prob = ODEProblem(SLC!, U0_, t_span, [t_0, k, r, K, H, d, Smax, gamma, cij])
        
        #Resolver el problema
        sol = solve(prob, maxiters=500)
        
        #Almacenar las soluciones de t, Na1, Na2, Sa1, Sa2
        for m in 1:size(sol.t , 1)
        push!(resultados_t, sol.t[m])
        push!(resultados_Na1, sol.u[m][1])
        push!(resultados_Na2, sol.u[m][2])
        push!(resultados_Sa1, sol.u[m][3])
        push!(resultados_Sa2, sol.u[m][4])
        end
        
        #Concatenar los resultados para obtener las soluciones completas de Na1, Na2, Sa1, Sa2
        #Iterations (days)
        resultados_t_concatenados = vcat(resultados_t...)
        #Abundance (individuals)
        resultados_Na1_concatenados = vcat(resultados_Na1...)
        resultados_Na2_concatenados = vcat(resultados_Na2...)
        #Size (mm)
        resultados_Sa1_concatenados = vcat(resultados_Sa1...)
        resultados_Sa2_concatenados = vcat(resultados_Sa2...)
    
        end  #fin del bucle de simulaciones

        #Graficar los resultados
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
    
plot!(legend=false)
xlims!(0, 7.0*10^4)
ylims!(0, 7.0*10^4)
plot!(background_color=:transparent, grid=true)


# #Fig 4b: with discrete leyend
# Variables para almacenar las frecuencias consolidadas
n_simulaciones = 100
n_bins = 100
consolidated_frequencies_ = zeros(Float64, n_bins, n_bins)
heat_0=surface()
j=1
#n=1

for j in 1
  cij = cij_span[j]  # Componente de competencia simétrica
  for n in 1:11
      H = H_span[n]  # Valor de explotación      
      # Reiniciar resultados para esta combinación de H y cij
      resultados_Na1_concatenados = Float64[]
      resultados_Na2_concatenados = Float64[]
      resultados_Sa1_concatenados = Float64[]
      resultados_Sa2_concatenados = Float64[]

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
              push!(resultados_Sa1_concatenados, sol.u[m][3])
              push!(resultados_Sa2_concatenados, sol.u[m][4])

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
        heat_0 = surface!(x_bins,
         y_bins, 
         consolidated_frequencies_,
         xlabel="Na1", 
         ylabel="Na2", 
         title="(Na1 vs Na2)", 
         color=cgrad(:thermal, rev=false),
         clims=(minimum(consolidated_frequencies_), maximum(consolidated_frequencies_)))
      else
        heat_0 = surface(x_bins,
         y_bins, 
         consolidated_frequencies_,
         xlabel="Na1", 
         ylabel="Na2", 
         title="(Na1 vs Na2)", 
         color=cgrad(:thermal, rev=false),
         clims=(minimum(consolidated_frequencies_), maximum(consolidated_frequencies_)))
      end
      display(heat_0)
  end
end

# Mostrar gráfico final
heat_0
surface!(background_color=:transparent, grid=true)
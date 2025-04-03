using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics
using Random
using Distributions
using StatsPlots
using DataFrames

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
    
    # Dynamics of the SLC
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
avg_oocytes = [77404, 385613] 
reggs = avg_oocytes ./ (365.14 * 0.42) 
r_ = reggs .* 0.998611 .* 0.971057 .* 0.4820525 .* 0.00629  

# Superposición reproductiva (mínimo esfuerzo reproductivo)
r_overlap = minimum(r_) / sum(r_) * sum(r_)

# Coeficientes de competencia
# c_12, c_21 = r_[1] / (sum(r_) - r_overlap), r_[2] / (sum(r_) - r_overlap)  #Interspecific  competition
c_11, c_22 = r_[1] / (r_[1] + r_overlap), r_[2] / (r_[2] + r_overlap)        #Intraspecific competition

# Parámetros biológicos
g_ = [0.998611, 0.971057, 0.4820525, 0.00629] 
d_ = [0.55, 0.59] ./ 365.14
size_growth_rate = [0.32, 0.36]
K_ = 64000.0 
k_ = 0.42
Smax_ = 53.0

# Intervalo y tasas de crecimiento
#N_span = 11
#H_span = range(0, 1, length=N_span) |> collect
#cij_span = range(0, 1, length=N_span) |> collect

N_span = 101
H_span = range(0, 1, length=N_span) |> collect
cij_span = range(0, 1, length=N_span) |> collect

gamma_ = [0.32, 0.36] ./ (365 * 0.42)

# Tamaños promedio de adultos antes y después de la implementacion de la normativa 
S_A_MPA_FA = [46.44, 44.45]
H_emp = [0.693, 0.57]

# -----------------------------
# Funciones auxiliares
# -----------------------------
function reproductive_period(t, t_0, k)
    return (t % t_0 / t_0) >= k ? 1.0 : 0.0
end

function reproductive_capacity(avg_size, Smax)
    Smat = 1.34 * avg_size - 28.06
    return min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax - Smat)), 0.0), 1.0)
end

# Modelos de escenarios
extinction_scenario() = (0.0, 0.0)

function one_species_extinction(cii,cjj, r, R, d, H, K, survivor)
    return survivor == 1 ? (K*(r[1] * R[1] - d[1] - H) / (cii + r[1] * R[1]), 0.0) :
                           (0.0, K*(r[2] * R[2] - d[2] - H) / (cjj + r[2] * R[2]))
end


function coexistence_scenario(cij, cji, cii, cjj, r, R, d, H, K)
if cij == 0.0 || cji == 0.0
    N1 = (K*(r[1] * R[1] - d[1] - H) / (r[1] * R[1]))
    N2 = (K*(r[2] * R[2] - d[2] - H) / (r[2] * R[2]))
    return (abs(N1), abs(N2))
else

Gamma = K^(-1)
rho1 = r[1] * R[1] - d[1] - H
rho2 = r[2] * R[2] - d[2] - H
z1 = cii + r[1] * R[1] * Gamma
z2 = cjj + r[2] * R[2] * Gamma
denom = cij * cji - z1 * z2

# Poblaciones de equilibrio
N1, N2 = (cij * rho1 - z1 * rho1) / denom, (cji * rho2 - z2 * rho2) / denom

# Asegurar que las poblaciones no sean negativas
return (max(N1), max(N2))
end
end


#Fig 4a: with discrete leyend
# -----------------------------
# Simulaciones
# -----------------------------
N_simulations = 10
t0_ = 365.14 * 0.42


extinction_results, one_species_results, coexistence_results = [], [], []


Threads.@threads for n in 1:N_span
    cij = cij_span[n]
    cji = cij
    for j in 1:N_span
    H = H_span[j]
    for i in 1:N_simulations
    k = k_ + k_*1/10 * randn()
    r = [r_[1] + (r_[1])*1/10 * randn()
    , r_[2] + (r_[2])*1/10 * randn()
    ] 
    K = K_ + (K_)*1/10 * randn()
    d = [d_[1] + (d_[1])*1/10* randn(), d_[2] + (d_[2])*1/10 * randn()]
    Smax = Smax_ + (Smax_)*1/10 * randn()
    
    R_ = [reproductive_capacity(S_A_MPA_FA[1], Smax), reproductive_capacity(S_A_MPA_FA[2], Smax)]
    
    push!(extinction_results, (cij, H, extinction_scenario()))
    push!(one_species_results, (cij, H, one_species_extinction(c_11,c_22, r, R_, d, H, K, 1), 
                                one_species_extinction(c_11,c_22, r, R_, d, H, K, 2)))
    push!(coexistence_results, (cij, H, coexistence_scenario(cij, cji,c_11,c_22, r_, R_, d, H, K)))
    end
end
end

extinction_results
one_species_results 
coexistence_results
# -----------------------------
# Resultados
# -----------------------------
df1 = DataFrame(cij= getindex.(extinction_results, 1),
                H = getindex.(extinction_results, 2),
                N_1 = getindex.(getindex.(extinction_results, 3),1),
                N_2 = getindex.(getindex.(extinction_results, 3),1))
df2 = DataFrame(cij = first.(one_species_results),
                H = getindex.(getindex.(one_species_results, 2), 1),
                N1_PRIMA = getindex.(getindex.(one_species_results, 3), 1),
                N2_0 = getindex.(getindex.(one_species_results, 3), 2),
                N1_0 = getindex.(getindex.(one_species_results, 4), 1),
                N2_PRIMA = getindex.(getindex.(one_species_results, 4), 2))
df3 = DataFrame(cij= getindex.(coexistence_results, 1),
                H = getindex.(coexistence_results, 2),
                N_1 = getindex.(getindex.(coexistence_results, 3),1),
                N_2 = getindex.(getindex.(coexistence_results, 3),2))
                show(df3, allrows=true)
# Filtrar valores positivos
df1_H_0_c0 = filter(row -> row.H == 0.0 && row.cij == 0.0, df1)
df1_H_0_c1 = filter(row -> row.H == 0.0 && row.cij == 1.0, df1)
df1_H_1_c0 = filter(row -> row.H == 1.0 && row.cij == 0.0, df1)
df1_H_1_c1 = filter(row -> row.H == 1.0 && row.cij == 1.0, df1)


df2_H_0_c0 = filter(row -> row.H == 0 && row.cij == 0, df2)
df2_H_0_c1 = filter(row -> row.H == 0 && row.cij == 1, df2)
df2_H_1_c0 = filter(row -> row.H == 1 && row.cij == 0, df2)
df2_H_1_c1 = filter(row -> row.H == 1 && row.cij == 1, df2)


df3_H_0_c0 = filter(row -> row.H == 0 && row.cij == 0, df3)
df3_H_0_c1 = filter(row -> row.H == 0 && row.cij == 1, df3)
df3_H_1_c0 = filter(row -> row.H == 1 && row.cij == 0, df3)
df3_H_1_c1 = filter(row -> row.H == 1 && row.cij == 1, df3)

df3_H_05_c0 = filter(row -> row.H == 0.5 && row.cij == 0, df3)
df3_H_0_c05 = filter(row -> row.H == 0 && row.cij == 0.5, df3)
df3_H_0_c05 = filter(row -> row.H == 0.5 && row.cij == 0.5, df3)

df3_H_0 = filter(row -> row.H == 0, df3)
df3_c_0 = filter(row -> row.cij == 0, df3)

df3_H_05 = filter(row -> row.H == 1, df3)
df3_c_055 = filter(row -> row.cij == 1, df3)


df3_H_1 = filter(row -> row.H == 1, df3)
df3_c_1 = filter(row -> row.cij == 1, df3)
 
# Coexistence scemarop
scatter(log10.(df3_H_0_c0.N_1), log10.(df3_H_0_c0.N_2), label="H = 0; cij = 0", color=:yellow, legend=:outertop)
scatter!(log10.(df3_H_0_c1.N_1), log10.(df3_H_0_c1.N_2), label="H = 0; cij = 1", color=:blue, legend=:outertop)
scatter!(log10.(df3_H_1_c0.N_1), log10.(df3_H_1_c0.N_2), label="H = 1; cij = 0", color=:green, legend=:outertop)
scatter!(log10.(df3_H_1_c1.N_1), log10.(df3_H_1_c1.N_2), label="H = 1; cij = 1", color=:red, legend=:outertop)

scatter!(log10.(df3_H_05_c0.N_1), log10.(df3_H_05_c0.N_2), label="H = 0.5; cij = 0", color=:black, legend=:outertop)
scatter!(log10.(df3_H_0_c05.N_1), log10.(df3_H_0_c05.N_2), label="H = 0; cij = 0.5", color=:brown, legend=:outertop)
scatter!(log10.(df3_H_05_c05.N_1), log10.(df3_H_05_c05.N_2), label="H = 0.5; cij = 0.5", color=:orange, legend=:outertop)

scatter!((df3_c_0.N_1), (df3_c_0.N_2), label="H=[0,1]; cij = 0", color=:orange, legend=:outertop)
scatter!((df3_H_0.N_1), (df3_H_0.N_2), label="cij= [0,1]; H = 0; ", color=:brown, legend=:outertop)
scatter!((df3_H_1.N_1), (df3_H_1.N_2), label="cij= [0,1]; H = 1; ", color=:red, legend=:outertop)
scatter((df3_c_1.N_1), (df3_c_1.N_2), label="H= [0,1]; cij = 1; ", color=:purple, legend=:outertop)


xlabel!("N1")
ylabel!("N2")
savefig("Fig4a.png")

n_bins = 101
mean_N1 = zeros(Float64, n_bins, n_bins)
mean_N2 = zeros(Float64, n_bins, n_bins)
# Filtrar los datos para H=0 y cij=0

# Iterar sobre todas las combinaciones de H y cij
for n in 1:N_span
    H = H_span[n]
    for j in 1:N_span
        cij = cij_span[j]
        
        # Filtrar los datos para la combinación actual de H y cij
        subset = filter(row -> row.H == H && row.cij == cij, df3)
        
        # Calcular la media de N_1 y N_2 para la combinación actual de H y cij
        if !isempty(subset)
            mean_N1[n, j] = mean(filter(x -> x > 0, subset.N_1))  # Ignorar valores no positivos
            mean_N2[n, j] = mean(filter(x -> x > 0, subset.N_2))  # Ignorar valores no positivos
        else
            mean_N1[n, j] = 0  # Asignar NaN si no hay datos
            mean_N2[n, j] = 0  # Asignar NaN si no hay datos
        end
    end
end
log10.(mean_N1) 
log10.(mean_N2) 

# Calcular los límites del rango de valores de N_1 y N_2 en el DataFrame original
min_N1 = log10(minimum(filter(x -> x > 0, df3.N_1))) # Ignorar valores no positivos
max_N1 = log10(maximum(df3.N_1))

min_N2 = log10(minimum(filter(x -> x > 0, df3.N_2)))  # Ignorar valores no positivos
max_N2 = log10(maximum(df3.N_2))





# Crear un gradiente de color continuo
color_gradient = cgrad(:viridis)  # Puedes cambiar a otras paletas como :plasma, :inferno, etc.

# Crear heatmaps con el gradiente basado en el rango de valores de df3
heatmap(H_span, cij_span, log10.(mean_N1),
        xlabel="H", ylabel="cij", title="N1",
        color=color_gradient, clims=(min_N1, max_N1))
savefig("N1_heatmap.png")
heatmap(H_span, cij_span, log10.(mean_N2),
        xlabel="H", ylabel="cij", title="N2",
        color=color_gradient, clims=(min_N2, max_N2))
savefig("N2_heatmap.png")


heatmap(df3.H, df3.cij, df3.N_1)

minimum(df)

# Filtrar valores positivos de N_1
df3_positive = filter(row -> row.N_1 > 0, df3)

# Obtener valores únicos de H y cij
unique_H = unique(df3_positive.H)
unique_cij = unique(df3_positive.cij)

# Crear una matriz para los valores de N_1
heatmap_matrix = fill(NaN, length(unique_H), length(unique_cij))

# Llenar la matriz con los valores de N_1
for row in eachrow(df3_positive)
    i = findfirst(==(row.H), unique_H)
    j = findfirst(==(row.cij), unique_cij)
    heatmap_matrix[i, j] = row.N_1
end

# Crear el heatmap
heatmap(unique_H, unique_cij, log10.(heatmap_matrix),
        xlabel="H", ylabel="cij", title="Heatmap de N1",
        color=:viridis, clims=(log10.(minimum(df3_positive.N_1)), log10.(maximum(df3_positive.N_1))))


       #Agrupar por las combinaciones únicas de cij y H
        grouped_df3 = groupby(df3, [:cij, :H])
        
        # Calcular el promedio de N_1 para cada combinación única de cij y H
        averaged_df3 = combine(grouped_df3, :N_1 => mean => :mean_N_1)
        
        # Obtener valores únicos de H y cij
        unique_H = unique(averaged_df3.H)
        unique_cij = unique(averaged_df3.cij)
        
        # Crear una matriz para los valores promedio de N_1
        heatmap_matrix = fill(NaN, length(unique_H), length(unique_cij))
        
        # Llenar la matriz con los valores promedio de N_1
        for row in eachrow(averaged_df3)
            i = findfirst(==(row.H), unique_H)
            j = findfirst(==(row.cij), unique_cij)
            heatmap_matrix[i, j] = row.mean_N_1
        end
        
        using NaNMath  # Para manejar NaN si es necesario

# Filtrar valores positivos en la matriz del heatmap
filtered_heatmap_matrix = map(x -> x > 0 ? x : NaN, heatmap_matrix)

# Crear el heatmap con la matriz filtrada
heatmap(unique_H, unique_cij, filtered_heatmap_matrix,
        xlabel="H", ylabel="cij", title="Heatmap de N1 (Solo Positivos)",
        color=:viridis, clims=(log10.(minimum(filtered_heatmap_matrix, dims=1)),
         log10.(maximum(filtered_heatmap_matrix, dims=1))))
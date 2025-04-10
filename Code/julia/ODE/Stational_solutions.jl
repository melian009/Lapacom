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
using NaNMath  


#=```
#Parameters and variables:
# - Ne = eggs abundance
# - Nt = trocophore abuncance
# - Nv = veliger abuncance
# - Nj = juvenile abuncance
# - Na = adults abuncance

# - Sa = adults size (Average sizes Before MPS+FULL = [46.44,44.45])
 
# - r = population growth rate [9.17,5.03]
# - R = reproductive capacity
# - K = carrying capacity (k = 1e^4)
# - X = Reproductive period [1,0] 
# - (1-X) = Exploitation period [1,0]
# - H = Exploitation rate (H = [0.639,0.57])
# - gEA = instant conversion rate of life stages (EA = Eggs to Adults) (gEA = 0.006)
# - de = natural mortality rate or death rate for eggs (de = [de_po,de_pa]) # Note: this value needs to be defined.
# - da = natural mortality rate or death rate for adults (da = [0.55,0.59]) # Note: empirical estimated values
# - Smax = maximum adult size estimated (56.0mm) # Note: Empirical value from the sample analized
# - gamma = adult growth rate (gamma=[0.32,0.36] year^{-1})
#```=#

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


  
#   Parameters for SLC 
avg_oocytes = [77404, 385613] 
reggs = avg_oocytes ./ (365.14 * 0.42) 
r_ = reggs .* 0.998611 .* 0.971057 .* 0.4820525 .* 0.00629  

# Superposición reproductiva (mínimo esfuerzo reproductivo)
r_overlap = minimum(r_) / sum(r_) * sum(r_)

# Coeficientes de competencia
c_12, c_21 = r_[1] / (sum(r_) - r_overlap), r_[2] / (sum(r_) - r_overlap)  #Interspecific  competition
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

# Filtrar los datos para H=0
df3_H_0 = subset(df3, :H => ByRow(==(0)))
df3_c_0 = subset(df3, :cij => ByRow(==(0)))


scatter(df3_c_0.N_1, df3_c_0.N_2,
        zcolor=df3_c_0[!, :H],  # Usar los valores de cij para el gradiente de color
        xlabel="N1", ylabel="N2", title="",
        colorbar_title="", palette=:plasma, label="cij=0, H = [0,1]")

scatter!(df3_H_0.N_1, df3_H_0.N_2,
        zcolor=df3_H_0[!, :cij],  # Usar los valores de cij para el gradiente de color
        xlabel="N1", ylabel="N2", title="",
        colorbar_title="H = [0,1]; cij = [0,1]", palette=:plasma,
        markerstrokecolor=:red,markerstrokewidth=1, label = "H=0, cij = [0,1]")

savefig("Figure_4a_grad.png")



n_bins = 101
mean_N1 = zeros(Float64, n_bins, n_bins)
mean_N2 = zeros(Float64, n_bins, n_bins)
# Filtrar los datos para H=0 y cij=0

# Iterar sobre todas las combinaciones de H y cij
# Iterar sobre todas las combinaciones de H y cij
for n in 1:N_span
    cij = cij_span[n]
    for j in 1:N_span
        H = H_span[j]

        # Filtrar los datos para la combinación actual de H y cij del DataFrame df3
        subset = filter(row -> hasproperty(row, :H) && hasproperty(row, :cij) && row.H == H && row.cij == cij, df3)

        # Verificar si el subset no está vacío
        if !isempty(subset)
            # Reemplazar valores negativos por ceros en las columnas N_1 y N_2
            subset.N_1 .= map(x -> x < 0 ? 0 : x, subset.N_1)
            subset.N_2 .= map(x -> x < 0 ? 0 : x, subset.N_2)

            # Calcular la media de los valores ajustados
            mean_N1[n, j] = mean(subset.N_1)
            mean_N2[n, j] = mean(subset.N_2)
        else
            # Si el subset está vacío, asignar ceros
            mean_N1[n, j] = 0.0
            mean_N2[n, j] = 0.0
        end
    end
end

mean_N1
mean_N2


# Calcular los límites del rango de valores de N_1 y N_2 en el DataFrame original
min_N1 = minimum(filter(x -> x > 0, df3.N_1)) # Ignorar valores no positivos
max_N1 = maximum(filter(x -> x < K_, df3.N_1))

min_N2 = minimum(filter(x -> x > 0, df3.N_2))  # Ignorar valores no positivos
max_N2 = maximum(filter(x -> x < K_, df3.N_2))



# Crear un gradiente de color continuo
color_gradient = cgrad(:viridis)

# Crear y guardar heatmaps para mean_N1 y mean_N2
for (data, title, filename, clims) in [(mean_N1, "N1", "N1_heatmap.png", (log10.(min_N1), log10.(max_N1))),
                                       (mean_N2, "N2", "N2_heatmap.png", (log10.(min_N2), log10.(max_N2)))]
    hm =heatmap(H_span, cij_span, log10.(data), xlabel="H", ylabel="cij", title=title, 
    color=color_gradient, clims=clims)
    display(hm)
end

# Agrupar por combinaciones únicas de cij y H y calcular promedios
averaged_df3_N1 = combine(groupby(df3, [:cij, :H]), :N_1 => mean => :mean_N_1)
averaged_df3_N2 = combine(groupby(df3, [:cij, :H]), :N_2 => mean => :mean_N_2)


# Crear matriz para los valores promedio de N_1
unique_H_N1, unique_cij_N1 = unique(averaged_df3_N1.H), unique(averaged_df3_N1.cij)
heatmap_matrix_N1 = fill(NaN, length(unique_cij_N1), length(unique_H_N1))
for row in eachrow(averaged_df3_N1)
    heatmap_matrix_N1[findfirst(==(row.cij), unique_cij_N1), findfirst(==(row.H),
    unique_H_N1)] = row.mean_N_1
end
# Crear matriz para los valores promedio de N_2
unique_H_N2, unique_cij_N2 = unique(averaged_df3_N2.H), unique(averaged_df3_N2.cij)
heatmap_matrix_N2 = fill(NaN, length(unique_cij_N2), length(unique_H_N2))
for row in eachrow(averaged_df3_N2)
    heatmap_matrix_N2[findfirst(==(row.cij), unique_cij_N2), findfirst(==(row.H),
    unique_H_N2)] = row.mean_N_2
end


# Filtrar valores positivos en la matriz del heatmap y reemplazar Inf por K_
filtered_heatmap_matrix_N1 = map(x -> isinf(x) ? K_ : (x > 0 ? x : 1), heatmap_matrix_N1)
filtered_heatmap_matrix_N2 = map(x -> isinf(x) ? K_ : (x > 0 ? x : 1), heatmap_matrix_N2)


# Crear y guardar el heatmap filtrado
using Plots
pyplot()  # Cambiar al backend PyPlot

heatmap(unique_H_N1, unique_cij_N1, log2.(filtered_heatmap_matrix_N1),
        xlabel="H", ylabel="cij", title="H vs Cij for N1", clabel="Ln(N1)",
        color=:viridis, clims=(minimum(log2.(filtered_heatmap_matrix_N1)),
         maximum(log2.(filtered_heatmap_matrix_N1))))
         
         savefig("Heatmap_N1.png")

heatmap(unique_H_N2, unique_cij_N2, log2.(filtered_heatmap_matrix_N2),
        xlabel="H", ylabel="cij", title="H vs Cij for N2 ",
        color=:viridis, clims=(minimum(log2.(filtered_heatmap_matrix_N2)),
         maximum(log2.(filtered_heatmap_matrix_N2))),clabel="Ln(N2)")
savefig("Heatmap_N2.png")


# Crear y guardar el heatmap filtrado
surface(unique_H_N1, unique_cij_N1, log2.(filtered_heatmap_matrix_N1),
        xlabel="H", ylabel="cij", title="H vs Cij for N1",
        color=:viridis, clims=(minimum(log2.(filtered_heatmap_matrix_N1)),
         maximum(log2.(filtered_heatmap_matrix_N1))),
         camera=(80, 10))
#savefig("Surface_N1.png")
surface(unique_H_N2, unique_cij_N2, log2.(filtered_heatmap_matrix_N2),
        xlabel="H", ylabel="cij", title="H vs Cij for N2 ",
        color=:viridis, clims=(minimum(log2.(filtered_heatmap_matrix_N2)),
        maximum(log2.(filtered_heatmap_matrix_N2))),
        camera=(80, 10))
#savefig("Surface_N2.png")


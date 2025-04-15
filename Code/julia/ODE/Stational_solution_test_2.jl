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

N_span = 101
H_span = range(0, 1, length=N_span) |> collect
cij_span = range(0, 1, length=N_span) |> collect

gamma_ = [0.32, 0.36] ./ (365 * 0.42)

# Tamaños promedio de adultos antes y después de la implementacion de la normativa 
S_A_MPA_FA = [46.44, 44.45]
H_emp = [0.693, 0.57]

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
#end

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


N_simulations = 1
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

df3 = DataFrame(cij= getindex.(coexistence_results, 1),
                H = getindex.(coexistence_results, 2),
                N_1 = getindex.(getindex.(coexistence_results, 3),1),
                N_2 = getindex.(getindex.(coexistence_results, 3),2)
                )
               
                
#n_bins = 101
#mean_N1 = zeros(Float64, n_bins, n_bins)
#mean_N2 = zeros(Float64, n_bins, n_bins)


#for n in 1:N_span
#    cij = cij_span[n]
#     for j in 1:N_span
#        H = H_span[j]

        # Filtrar los datos para la combinación actual de H y cij del DataFrame df3
#        subset = filter(row -> hasproperty(row, :H) && hasproperty(row, :cij) && row.H == H && row.cij == cij, df3)

        # Verificar si el subset no está vacío
#        if !isempty(subset)
            # Reemplazar valores negativos por ceros en las columnas N_1 y N_2
#            subset.N_1 .= map(x -> x < 0 ? 0 : x, subset.N_1)
#            subset.N_2 .= map(x -> x < 0 ? 0 : x, subset.N_2)

            # Calcular la media de los valores ajustados
#            mean_N1[n, j] = mean(subset.N_1)
#            mean_N2[n, j] = mean(subset.N_2)
#        else
            # Si el subset está vacío, asignar ceros
#            mean_N1[n, j] = 0.0
#            mean_N2[n, j] = 0.0
#        end
#    end
#end


#size(df3.N_1)
#mean_N1
#mean_N2



# Calcular los límites del rango de valores de N_1 y N_2 en el DataFrame original
#min_N1 = minimum(filter(x -> x > 0, df3.N_1)) # Ignorar valores no positivos
#max_N1 = maximum(filter(x -> x < K_, df3.N_1))

#min_N2 = minimum(filter(x -> x > 0, df3.N_2))  # Ignorar valores no positivos
#max_N2 = maximum(filter(x -> x < K_, df3.N_2))


sort!(df3, [:cij, :H])
cij_vals = unique(df3.cij)
H_vals = unique(df3.H)
# Reshape into matrix N1
z = reshape(df3.N_1, length(H_vals), length(cij_vals))  # rows = H, cols = cij
heatmap(cij_vals,H_vals,z,xlabel="cij",ylabel="H",title="Heatmap of N_1",colorbar_title="N_1",c=:viridis)


# Reshape into matrix N2
#z = reshape(df3.N_2, length(H_vals), length(cij_vals))  # rows = H, cols = cij
#heatmap(cij_vals,H_vals,z,xlabel="cij",ylabel="H",title="Heatmap of N_2",colorbar_title="N_2",c=:viridis)


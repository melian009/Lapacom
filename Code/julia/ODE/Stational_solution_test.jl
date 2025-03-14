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

# -----------------------------
# Parámetros iniciales
# -----------------------------
avg_oocytes = [77404, 385613] 
reggs = avg_oocytes ./ (365.14 * 0.42) 
r_ = reggs .* 0.998611 .* 0.971057 .* 0.4820525 .* 0.00629  

# Superposición reproductiva (mínimo esfuerzo reproductivo)
r_overlap = minimum(r_) / sum(r_) * sum(r_)

# Coeficientes de competencia
# c_12, c_21 = r_[1] / (sum(r_) - r_overlap), r_[2] / (sum(r_) - r_overlap)
c_11, c_22 = r_[1] / (r_[1] + r_overlap), r_[2] / (r_[2] + r_overlap)

# Parámetros biológicos
g_ = [0.998611, 0.971057, 0.4820525, 0.00629] 
d_ = [0.55, 0.59] ./ 365.14
size_growth_rate = [0.32, 0.36]
K_ = 64000.0 
k_ = 0.42
Smax_ = 53.0

# Intervalo y tasas de crecimiento
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
    k, r, K = k_ + 0.1 * randn(), [r_[1] + 0.1 * randn(), r_[2] + 0.1 * randn()], K_ + 0.1 * randn()
    Smax = Smax_ + randn()
    R_ = [reproductive_capacity(S_A_MPA_FA[1], Smax), reproductive_capacity(S_A_MPA_FA[2], Smax)]
    
    push!(extinction_results, (cij, H, extinction_scenario()))
    push!(one_species_results, (cij, H, one_species_extinction(c_11,c_22, r, R_, d_, H, K_, 1), 
                                one_species_extinction(c_11,c_22, r, R_, d_, H, K_, 2)))
    push!(coexistence_results, (cij, H, coexistence_scenario(cij, cji,c_11,c_22, r_, R_, d_, H, K_)))
    end
end

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

# Filtrar valores positivos
df1_positive = filter(row -> row.N_1 >= 0 && row.N_2 >= 0, df1)
df2_positive = filter(row -> row.N1_PRIMA > 1 && row.N2_PRIMA > 1, df2)
df3_positive = filter(row -> row.N_1 > 0 && row.N_2 > 0, df3)
df3_NO_0 = filter(row -> row.H > 0 && row.cij > 0, df3)

# -----------------------------
# Visualización
# -----------------------------
scatter(df2.N1_PRIMA, df2.N2_0, label="N2 = 0;  N1 = N1*", color=:red)
scatter!(df2.N1_0, df2.N2_PRIMA, label="N1 = 0;  N2 = N2*", color=:green)
scatter!(df3.N_1, df3.N_2, label="N1 = N1*; N2 = N2*", color=:blue, legend=:outertop)
scatter!(df1.N_1, df1.N_2, label="N1 = N2 = 0", color=:black)
xlabel!("N1")
ylabel!("N2")


# Exploitation Gradient
scatter(df1_positive.H, df1_positive.N_1, label="N1 = 0", color=:black)
scatter!(df1_positive.H, df1_positive.N_2, label="N2 = 0", color=:grey)
scatter!(df2_positive.H, df2_positive.N1_PRIMA, label="N2 = 0;  N1 = N1*", color=:red)
scatter!(df2_positive.H, df2_positive.N2_PRIMA, label="N1 = 0;  N2 = N2*", color=:pink)
scatter!(df3_positive.H, df3_positive.N_1, label="N1 = N1* Coexistencia", color=:darkblue)
scatter!(df3_positive.H, df3_positive.N_2, label="N2 = N2* Coexistencia", color=:blue,legend=:outertop)
ylims!(0,100)

# Competence gradient
scatter(df1.cij, df1.N_1, label="N1 = 0", color=:black)
scatter!(df1.cij, df1.N_2, label="N2 = 0", color=:grey)
scatter!(df2.cij, df2.N1_PRIMA, label="N2 = 0;  N1 = N1*", color=:red)
scatter!(df2.cij, df2.N2_PRIMA, label="N1 = 0;  N2 = N2*", color=:pink)
scatter!(df3.cij, df3.N_1, label="N1 = N1* Coexistencia", color=:darkblue)
scatter!(df3.cij, df3.N_2, label="N2 = N2* Coexistencia", color=:yellow)
ylims!(0,100)


# Discrete colors for competence and exploitation range values. 
# Mapping different colors for unice values of cij and H.
cij_values = unique(df1.cij)
H_values = unique(df1.H)

# Mapping each unique combination of cij and H to a distinct color.
color_map = Dict()

# Function to assign a color to each combination (cij, H)
function get_color(cij, H)
    key = (cij, H)
    if !haskey(color_map, key)
        #If it doesn't have a color assigned, it generates a random one.
        color_map[key] = RGB(rand(), rand(), rand())
    end
    return color_map[key]
end

# Function to obtain opacity (based on cij and H)
function get_opacity(cij, H)
    # Normalize the cij and H values ​​to obtain an opacity between 0 and 1.
    cij_norm = (cij - minimum(cij_values)) / (maximum(cij_values) - minimum(cij_values))
    H_norm = (H - minimum(H_values)) / (maximum(H_values) - minimum(H_values))
    # Use an average of the normalization for opacity.
    return (cij_norm + H_norm) / 2
end


# Plot the points of the first scenario (df1_positive) with colors and opacity according to cij and H
scatter(df1.N_1, df1.N_2, label="N1 = N2 = 0", 
        xlabel="N1", ylabel="N2",
        markers=(:diamond, 5),
        color=[get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df1.cij, df1.H)])

# Only one survivor - Patella ordinaris (df2_positive)
scatter!(df2.N1_PRIMA, df2.N2_0, 
         label="N2 = 0;  N1 = N1*", 
         color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        markers=(:circle, 5))


# Only one survivor - Patella aspera (df2_positive)
scatter!(df2.N1_0, df2.N2_PRIMA, 
         label="N1 = 0;  N2 = N2*", 
         color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         markers=(:square, 5))

# Coexistence (df3_positive)
scatter!(df3.N_1, df3.N_2, 
         label="N1 = N1* ;  N2 = N2*", 
         color=[get_color(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
markers=(:hexagon, 5), legend=:outertop)





# Exploitation Gradient
scatter(df1.H, df1.N_1, label="N1 = 0", 
color=[get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        markers=(:utriangle, 5))

scatter!(df1.H, df1.N_2, label="N2 = 0", 
color=[get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        markers=(:dtriangle, 5))

scatter!(df2.H, df2.N1_PRIMA, label="N2 = 0;  N1 = N1*", 
color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        markers=(:square, 5))


scatter!(df2.H, df2.N2_PRIMA, label="N1 = 0;  N2 = N2*",
        color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        markers=(:hentagon, 5))


scatter!(df3.H, df3.N_1, label="N1 = N1* Coexistencia",
 color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
 marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
 markers=(:diamond, 5))


scatter(df3.H, df3.N_2, label="N2 = N2* Coexistencia",  color=[get_color(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
 marker_z=[get_opacity(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
markers=(:pentagon, 5), legend=:outertop)

ylims!(0,40)

# Competence gradient
scatter(df1.cij, df1.N_1, label="N1 = 0", 
color=[get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        markers=(:utriangle, 5))

scatter!(df1.cij, df1.N_2, label="N2 = 0", 
color=[get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        markers=(:dtriangle, 5))

scatter!(df2.cij, df2.N1_PRIMA, label="N2 = 0;  N1 = N1*", 
color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        markers=(:square, 5))


scatter!(df2.cij, df2.N2_PRIMA, label="N1 = 0;  N2 = N2*",
        color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        markers=(:hentagon, 5))


scatter!(df3.cij, df3.N_1, label="N1 = N1* Coexistencia",
 color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
 marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
 markers=(:diamond, 5))


scatter!(df3.cij, df3.N_2, label="N2 = N2* Coexistencia",  color=[get_color(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
 marker_z=[get_opacity(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
markers=(:pentagon, 5), legend=:outertop)

ylims!(0,10)



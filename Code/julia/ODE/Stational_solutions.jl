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
 - Ne = eggs abundance (CLC)
 - Nt = trocophore abuncance (CLC)
 - Nv = veliger abuncance (CLC)
 - Nj = juvenile abuncance (CLC)
 - Na = adults abuncance (SLC and CLC)

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

 SLC simplified ecuation
 dN_i/dt = r_i * R_i * N_i * ((K - N_i)/ K) - d_i * N_i - (1 - periodX(t)) * H * N_i - N_i * sum(c_ij * N_j)
 
 Stability Scenarios for:
 N1 = Patella ordinaria
 N2 = Patella aspera

 Total extiontion scenarios
 N1* = N2* = 0

 One species extinction scenario
 Patella ordinaria survives
 N2* = 0; N1* = (r_1*R_1 - d_1 - H_1*(1-X))/(r_1*R_1*gamma)
 
 Patella aspera survives
 N1* = 0; N2* = (r_2*R_2 - d_2 - H_2*(1-X))/(r_2*R_2*gamma)

 Where:
 z_i = c_ij + r_i*R_i*gamma
 rho_i = r_i*R_i - d_i - H*(1-X)
 gamma = K^{-1}

 Cohexistance scenario cij = 0
 N1* = (c21*rho1 - z1*rho2)/(c12*c21 - z1*z2)
 N2* = (c12*rho2 - z2*rho1)/(21*c12 - z2*z1)
```

# Parámetros iniciales
avg_oocytes = [77404, 385613]
reggs = avg_oocytes ./ (365 * 0.42)
r_ = reggs .* 0.998611 .* 0.971057 .* 0.4820525 .* 0.00629
d_ = [0.55, 0.59] ./ 365.14
size_growth_rate = [0.32, 0.36]
k_ = 0.42
K_ = 64000.0
Smax_ = 53.0
N_span = 11  # Número de valores
H_r = range(0, 1, length=N_span)
Cij_r = range(0, 1, length=N_span)
H_span = collect(H_r)
cij_span = collect(Cij_r)



# Funciones para calcular soluciones numéricas de los escenarios de estabilidad
function extinction_scenario()
    return (0.0, 0.0)
end

function one_species_extinction_scenario(r, R, d, H, Gamma, survivor)
    if survivor == 1
        return ((r[1] * R[1] - d[1] - H) / (r[1] * R[1] * Gamma), 0.0)
    else
        return (0.0, (r[2] * R[2] - d[2] - H) / (r[2] * R[2] * Gamma))
    end
end

function coexistence_scenario(cij, r, R, d, H, Gamma)
    if cij == 0
        N1 = ((r[1] * R[1] - d[1] - H) / (r[1] * R[1] * Gamma), 0.0)
        N2 = (0.0, (r[2] * R[2] - d[2] - H) / (r[2] * R[2] * Gamma))
        return (N1, N2)
    else
    rho1 = r[1] * R[1] - d[1] - H
    rho2 = r[2] * R[2] - d[2] - H 
    z1 = cij + r[1] * R[1] * Gamma
    z2 = cij + r[2] * R[2] * Gamma
    
    N1 = (cij * rho1 - z1 * rho2) / (cij*cji - z1 * z2)
    N2 = (cji * rho2 - z2 * rho1) / (cji*cij - z1 * z2)
    return (N1, N2)
    end
end

# Simulación de escenarios
extinction_results = []
one_species_results = []
coexistence_results = []

N_simulations = 500

# Bucle que recorre solo un valor de H y un valor de cij por vuelta
for j in 1:N_span  # Esto es para recorrer los valores de cij
    cij = cij_span[j]
    
    for h in 1:N_span  # Esto es para recorrer los valores de H
        H = H_span[h]
        for i in 1:N_simulations
        # Parámetros con ruido
        k = k_ + 0.1 * randn()
        r = [r_[1] + 0.1 * randn(), r_[2] + 0.1 * randn()] 
        K = K_ + 0.1 * randn()
        gamma = [size_growth_rate[1] + 0.1 * randn(), size_growth_rate[2] + 0.1 * randn()]
        d = [d_[1], d_[2]]
        Smax = Smax_ + 1 * randn()
        Gamma = 1 / K_
        
        # Obtener soluciones de cada escenario
        ext = extinction_scenario()
        patella_ord_survives = one_species_extinction_scenario(r_, [0.5, 0.5], d_, H, Gamma, 1)
        patella_asp_survives = one_species_extinction_scenario(r_, [0.5, 0.5], d_, H, Gamma, 2)
        coexist = coexistence_scenario(cij, r_, [0.5, 0.5], d_, H, X, Gamma)
        
        # Almacenar los resultados para N1 y N2 en cada escenario
        push!(extinction_results, (cij, H, ext))
        push!(one_species_results, (cij, H, patella_ord_survives, patella_asp_survives))
        push!(coexistence_results, (cij, H, coexist))
        end
    end
end

# Visualizar resultados en tablas
 # extinction scenario
 df1 = DataFrame(cij = [r[1] for r in extinction_results], 
                H = [r[2] for r in extinction_results], 
                N_1 = [r[3][1] for r in extinction_results], 
                N_2 = [r[3][2] for r in extinction_results]) 
 show(df1, allrows=true) 
 # Filter positive values
 df1_positive = df1[df1.N_1 .>= 0 .&& df1.N_2 .>= 0, :]
 show(df1_positive, allrows=true) 

 # One species extinction scenario
df2 = DataFrame(cij = [r[1] for r in one_species_results], 
                H = [r[2] for r in one_species_results], 
                N1_PRIMA = [r[3][1] for r in one_species_results], 
                N2_0 = [r[3][2] for r in one_species_results],
                N1_0 = [r[4][1] for r in one_species_results], 
                N2_PRIMA = [r[4][2] for r in one_species_results])
 show(df2, allrows=true)
 #Filter positive values
 df2_positive = df2[df2.N1_PRIMA .>= 0 .&& df2.N2_0 .>= 0 .&& df2.N1_0 .>= 0 .&& df2.N2_PRIMA .>= 0, :]
 show(df2_positive, allrows=true)

# Coexistence scenario
df3 = DataFrame(cij = [r[1] for r in coexistence_results], 
                H = [r[2] for r in coexistence_results], 
                N_1 = [r[3][1] for r in coexistence_results], 
                N_2 = [r[3][2] for r in coexistence_results])
 show(df3, allrows=true)
 # Filter positive values
 df3_positive = df3[df3.N_1 .>= 0 .&& df3.N_2 .>= 0, :]
 show(df3_positive, allrows=true)


#Graficar los valores filtrados

# Extinción
scatter(df1_positive.cij, df1_positive.H, label="Extinción", xlabel="N1", ylabel="N2", color=:black)

# Un solo sobreviviente - Patella ordinaria
scatter!([r.N1_PRIMA for r in eachrow(df2_positive)], [r.N2_0 for r in eachrow(df2_positive)], label="N1 Prima (Patella ordinaria)", color=:blue)

# Un solo sobreviviente - Patella aspera
scatter!([r.N1_0 for r in eachrow(df2_positive)], [r.N2_PRIMA for r in eachrow(df2_positive)], label="N1 (Patella aspera)", color=:red)

# Coexistencia
scatter!([r.N_1 for r in eachrow(df3_positive)], [r.N_2 for r in eachrow(df3_positive)], label="Coexistencia", color=:green, legend=:outerright)



# Colores discretos para diferentes valores de cij y H. Aquí mapeamos los valores únicos de cij y H a colores distintos.
cij_values = unique(df1.cij)
H_values = unique(df1.H)

# Mapeamos cada combinación única de cij y H a un color distinto.
color_map = Dict()

# Función para asignar un color a cada combinación (cij, H)
function get_color(cij, H)
    key = (cij, H)
    if !haskey(color_map, key)
        # Si no tiene un color asignado, generamos uno aleatorio
        color_map[key] = RGB(rand(), rand(), rand())
    end
    return color_map[key]
end

# Función para obtener la opacidad (en base a cij y H)
function get_opacity(cij, H)
    # Normalizamos los valores de cij y H para obtener una opacidad entre 0 y 1.
    cij_norm = (cij - minimum(cij_values)) / (maximum(cij_values) - minimum(cij_values))
    H_norm = (H - minimum(H_values)) / (maximum(H_values) - minimum(H_values))
    # Usamos un promedio de la normalización para la opacidad.
    return (cij_norm + H_norm) / 2
end


# Graficar los puntos del primer escenario (df1_positive) con colores y opacidad según cij y H
scatter(df1.cij, df1.H, label="Total extinction", 
        xlabel="N1", ylabel="N2",
        markers=(:diamond, 5),
        color=[get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df1.cij, df1.H)])

# Un solo sobreviviente - Patella ordinaria (df2_positive)
scatter!([r.N1_PRIMA for r in eachrow(df2)], 
         [r.N2_0 for r in eachrow(df2)], 
         label="N1* (Patella ordinaria)", 
         color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        markers=(:circle, 5))


# Un solo sobreviviente - Patella aspera (df2_positive)
scatter!([r.N1_0 for r in eachrow(df2)], 
         [r.N2_PRIMA for r in eachrow(df2)], 
         label="N1 (Patella aspera)", 
         color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         markers=(:square, 5))

# Coexistencia (df3_positive)
scatter!([r.N_1 for r in eachrow(df3)], 
         [r.N_2 for r in eachrow(df3)], 
         label="Coexistencia", 
         color=[get_color(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
markers=(:hexagon, 5))

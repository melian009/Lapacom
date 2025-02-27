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

# Simple Life Cycle (SLC) model

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


# Parámetros iniciales [Na1, Na2]
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


# Stability scenarios

# N1 = N2 = 0
function extinction_scenario() 
    return (0.0, 0.0)
end

# N1 = 0; N2 = N2* or N2 = 0; N1 = N1* 
function one_species_extinction_scenario(cij, cji, r, R, d, H, Gamma, survivor)
    if survivor == 1
        return ((r[1] * R[1] - d[1] - H) / (cij + r[1] * R[1] * Gamma), 0.0)
    else
        return (0.0, (r[2] * R[2] - d[2] - H) / (cji + r[2] * R[2] * Gamma))
    end
end

#N1 = N1*; N2 = N2* || Cij = 0 & Cij =! 0
function coexistence_scenario(cij, cji, r, R, d, H, Gamma)
    if cij == 0 # Cij = 0
        N1 = ((r[1] * R[1] - d[1] - H) / (r[1] * R[1] * Gamma))
        N2 = ((r[2] * R[2] - d[2] - H) / (r[2] * R[2] * Gamma))
        return (N1, N2)

    else # Cij != 0
    rho1 = r[1] * R[1] - d[1] - H
    rho2 = r[2] * R[2] - d[2] - H 
    z1 = cij + r[1] * R[1] * Gamma
    z2 = cij + r[2] * R[2] * Gamma
    N1 = (cij * rho1 - z1 * rho2) / (cij*cji - z1 * z2)
    N2 = (cji * rho2 - z2 * rho1) / (cji*cij - z1 * z2)
    return (N1, N2)
    end
end

# Simulated scenarios
extinction_results = []
one_species_results = []
coexistence_results = []

N_simulations = 500
t_span = (0.0, 365.14*10)  # Tiempo de simulación (por ejemplo, un año)
t_plt = 0.0:1.0:365.14*10  # Los tiempos en los que se evaluará la solución
t0_ = 365.14*0.42

# Loops for descrete values of cij and H
for j in 1:N_span  # Run cij values
    cij = cij_span[j]
    cji = (cij)
    for h in 1  # Run H values
        H = H_span[h]
        #r i in 1:N_simulations
        # Noised parameters
        k = k_ + 0.1 * randn()
        r = [r_[1] + 0.1 * randn(), r_[2] + 0.1 * randn()] 
        K = K_ + 0.1 * randn()
        gamma = [size_growth_rate[1] + 0.1 * randn(), size_growth_rate[2] + 0.1 * randn()]
        d = [d_[1], d_[2]]
        Smax = Smax_ + 1 * randn()
        Gamma = 1 / K_
        

        # N_1
        avg_size_1 = 46.44 
        Smat_1 = 1.34 * (avg_size_1) - 28.06
        R_1 = min(max(0.5 * (1.0 + (avg_size_1 - Smat_1) / (Smax - Smat_1)), 0.0), 1.0)
        
        # N_2
        avg_size_2 = 44.45
        Smat_2 = 1.34 * (avg_size_1) - 28.06
        R_2 = min(max(0.5 * (1.0 + (avg_size_2 - Smat_2) / (Smax - Smat_2)), 0.0), 1.0)
        R_ = [R_1, R_2]
        # Solutions for each scenario
        ext = extinction_scenario()
        patella_ord_survives = one_species_extinction_scenario(cij, cji, r_, R_, d_, H, Gamma, 1)
        patella_asp_survives = one_species_extinction_scenario(cij, cji, r_, R_, d_, H, Gamma, 2)
        coexist = coexistence_scenario(cij, cji, r_, [0.5, 0.5], d_, H, Gamma)
        
        # Store the results for N1 and N2 in each scenario
        push!(extinction_results, (cij, H, ext))
        push!(one_species_results, (cij, H, patella_ord_survives, patella_asp_survives))
        push!(coexistence_results, (cij, H, coexist))
        #end
    end
end

# Display results in tables
 # Extinction scenario (N1 = N2 = 0)
 df1 = DataFrame(cij = [r[1] for r in extinction_results], 
                H = [r[2] for r in extinction_results], 
                N_1 = [r[3][1] for r in extinction_results], 
                N_2 = [r[3][2] for r in extinction_results]) 
 show(df1, allrows=true)

 # Filter positive values 
 df1_positive = df1[df1.N_1 .>= 0 .&& df1.N_2 .>= 0, :]
 show(df1_positive, allrows=true) 

 # One species extinction scenario (N2 = 0; N1 = N1* or N1 = 0; N2 = N2*)
 df2 = DataFrame(cij = [r[1] for r in one_species_results], 
                H = [r[2] for r in one_species_results], 
                N1_PRIMA = [r[3][1] for r in one_species_results], 
                N2_0 = [r[3][2] for r in one_species_results],
                N1_0 = [r[4][1] for r in one_species_results], 
                N2_PRIMA = [r[4][2] for r in one_species_results])
 
 #Filter positive values
 df2_positive = df2[df2.N1_PRIMA .>= 0 .&& df2.N2_0 .>= 0 .&& df2.N1_0 .>= 0 .&& df2.N2_PRIMA .>= 0, :]
 show(df2_positive, allrows=true)

# Coexistence scenario (N1 = N1* and N2 = N2*)
df3 = DataFrame(cij = [r[1] for r in coexistence_results], 
                H = [r[2] for r in coexistence_results], 
                N_1 = [r[3][1] for r in coexistence_results], 
                N_2 = [r[3][2] for r in coexistence_results])
 
# Filter positive values
 df3_positive = df3[df3.N_1 .>= 0 .&& df3.N_2 .>= 0, :]
 #show(df3, allrows=true)

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



#= Limit Cycles without leyend
limt_cycle = plot()
for j in 1:10:11 # Cij
    cij = cij_span[j]  #Simetric competence component
    for n in 1:11 # H 
        H = H_span[n] #Exploitation rate

        # Storing the simulation result sets
        resultados_t_concatenados = Float64[]  # To store the values ​​of t
        
        resultados_Na1_concatenados = Float64[]  # To store the values ​​of Na1
        resultados_Na2_concatenados = Float64[]  # To store the values ​​of Na2
        resultados_Sa1_concatenados = Float64[]  # To store the values ​​of Sa1
        resultados_Sa2_concatenados = Float64[]  # To store the values ​​of Sa2
    
        # Storing the simulation result sets
        resultados_t = Float64[]  # To store the values ​​of t
        
        resultados_Na1 = Float64[]  # To store the values ​​of Na1
        resultados_Na2 = Float64[]  # To store the values ​​of Na2
        resultados_Sa1 = Float64[]  # To store the values ​​of Sa1
        resultados_Sa2 = Float64[]  # To store the values ​​of Sa2
    
        # Simulations
        for i in 1:N_simulations
          # Generate random values ​​for parameters from a normalized distribution
          t_0 = t0_ + 0.0001 * randn()
          k = k_ + 0.01 * randn()
          r = [r_[1] + 0.01 * randn(), r_[2] + 0.01 * randn()] 
          K = K_ + 0.1 * randn()
          gamma = [size_growth_rate[1] + 0.01 * randn(),size_growth_rate[2] + 0.01 * randn()]
          d = [d_[1], d_[2]]
          Smax = Smax_ + 0.1 * randn()
        
          # Initial conditions
          U0_ = [10^4,10^4, mean([33.4,37.4]), mean([34.6,37.5])]
        
          # Define ODE problem
          prob = ODEProblem(SLC!, U0_, t_span, [t_0, k, r, K, H, d, Smax, gamma, cij])
        
          # Solve ODE problem
          sol = solve(prob, maxiters=1000)
        
          # To store solutions of t, Na1, Na2, Sa1, Sa2
          for m in 1:size(sol.t , 1)
            push!(resultados_t, sol.t[m])
            push!(resultados_Na1, sol.u[m][1])
            push!(resultados_Na2, sol.u[m][2])
            push!(resultados_Sa1, sol.u[m][3])
            push!(resultados_Sa2, sol.u[m][4])
          end
        
          # Concatenate the results to obtain the complete solutions of Na1, Na2, Sa1, Sa2
    
             resultados_t_concatenados = vcat(resultados_t...)
        
             resultados_Na1_concatenados = vcat(resultados_Na1...)
             resultados_Na2_concatenados = vcat(resultados_Na2...)
             resultados_Sa1_concatenados = vcat(resultados_Sa1...)
             resultados_Sa2_concatenados = vcat(resultados_Sa2...)
    
        end 
        # Plot results 
        if j == 1 && n == 1
          limt_cycle = plot!(resultados_Na1_concatenados, resultados_Na2_concatenados, 
                 xlabel="N1", ylabel="N2", 
                 label="H=$H, Cij=$cij",  color = [get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)])
        else
          limt_cycle =  plot!(resultados_Na1_concatenados, resultados_Na2_concatenados, 
                  label="H=$H, Cij=$cij", color = [get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)])
        end
    end
end
display(limt_cycle) =#

# Plot stability scenarios over the limit cycles
# Plot the points of the first scenario (df1_positive) with colors and opacity according to cij and H
scatter(df1_positive.N_1, df1_positive.N_2, label="N1 = N2 = 0", 
        xlabel="N1", ylabel="N2",
        markers=(:diamond, 5),
        color=[get_color(cij, H) for (cij, H) in zip(df1.cij, df1.H)],
        marker_z=[get_opacity(cij, H) for (cij, H) in zip(df1.cij, df1.H)])

# Only one survivor - Patella ordinaris (df2_positive)
scatter!([r.N1_PRIMA for r in eachrow(df2_positive)], 
         [r.N2_0 for r in eachrow(df2_positive)], 
         label="N2 = 0;  N1 = N1*", 
         color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
        markers=(:circle, 5))


# Only one survivor - Patella aspera (df2_positive)
scatter!([r.N1_0 for r in eachrow(df2_positive)], 
         [r.N2_PRIMA for r in eachrow(df2_positive)], 
         label="N1 = 0;  N2 = N2*", 
         color=[get_color(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df2.cij, df2.H)],
         markers=(:square, 5))

# Coexistence (df3_positive)
scatter!([r.N_1 for r in eachrow(df3_positive)], 
         [r.N_2 for r in eachrow(df3_positive)], 
         label="N1 = N1* ;  N2 = N2*", 
         color=[get_color(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
         marker_z=[get_opacity(cij, H) for (cij, H) in zip(df3.cij, df3.H)],
markers=(:hexagon, 5), legend=:outerright)


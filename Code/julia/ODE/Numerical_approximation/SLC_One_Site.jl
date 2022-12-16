## Packeges
using Pkg
Pkg.activate(".")
using LinearAlgebra
# using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
#using CairoMakie
using Statistics
using DataFrames
# using CSV
using DiffEqParamEstim
using Optim
using Plots
using Plots.PlotMeasures



### ------------------------------------------------------------------------------------------------
### 1. Model formulation for Simplex Life Cycle (SLC) : 3 model formulation
### ------------------------------------------------------------------------------------------------
"""
Variables:
  - Nⱼ: Egg Abundances
  - Nₐ: Adult Abundance
  - Sₐ: Adult Size
Parameters:
  - r: intrinsic growth rate
  - g: instant conversion rate of eggs turning into adult
  - d: death rate
  - E: exploitation rate
  - K: carrying capacity
  - X: Reproductive Cycle
"""

#= 
In the first formulation of the equations abundances aren't linked with the size.
=#

# The model formulation is:

function single_site!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = p
  du[1] = dNⱼ = (r * Nₐ * ((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
end

#=
The last equation doesn't consider a cycle of exploitation/reproduction. This supose that
the exploitation occurs continuously along Time
=#

# -----------------------------------------------------------------------------------------

#=
In the second definition equations are defined in such a way that the number of eggs (N_e) 
depends on the size of the adults (S_a) and the adult abundance (N_a) which depends on the
exploitation rate (Et) but still to not consider a reproductive/exploitation
cycle. This means that the three ODE are linked between each other but with a continuous 
expoitation with time.
=#

# The model formulation is:

function single_site_S!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ = p
  du[1] = dNⱼ = (r * Nₐ * (Sₐ/sizeₘₐₓ) * ((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (E(t) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ - (sizeₘₐₓ * E(t))))
end

#=
The first equation present the size component that link abundances equiations with size 
ODE. Also the last equation doesn't consider a cycle of exploitation/reproduction. This
supose that the exploitation occurs continuously alongt time.
=#

# -----------------------------------------------------------------------------------------

#= 
In the third definition equiations consider the reproductive/exploitation cycle (X) 
that change during the year (0 for explitation period and 1 for resproduction period).
=#

# The model definition is:

function single_site_S_X!(du, u, p, t)
  Nⱼ, Nₐ, Sₐ = u
  r, g, dⱼ, dₐ, Exp, K, size_growth_rate, sizeₘₐₓ, X, rate = p
  du[1] = dNⱼ = (X(t) * r * Nₐ * (Sₐ/sizeₘₐₓ)*((K - Nₐ) / K)) - (dⱼ * Nⱼ) - (g * Nⱼ)
  du[2] = dNₐ = (g * Nⱼ) - (dₐ * Nₐ) - (Exp(t,rate) * Nₐ)
  du[3] = dSₐ = size_growth_rate * Sₐ * (1 - Sₐ / (sizeₘₐₓ * (1 - Exp(t,rate))))
end
#=
In this formulation eggs abundances and adult size now have the component of the 
reproductiveX)/exploitation(1-X) cicle. This condition de model to give egg 
inputs when theres no exploitation, and adutls are exploited when are not producing eggs. 
=#

# The time varing exploitation is:

Et(t) = (sin(t)^2) / 2  

#=
To adjust the sin function horizonally (stretch/shrink the wave length), multiply t by a factor c.
You may also use a function with if statements:
=#

function Et(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return 0.5          # Proportion of the adult population exploited: 50%
    
  end
end 



### -----------------------------------------------------------------------------------
### 2. Graphic representation of SLC in a single site for the three model formulations.
### -----------------------------------------------------------------------------------

"""
Fist definition: 
Abundances ODE and Size ODE are not linked and no reproductive cicle consideration:
- With function SLC_single_site!
"""

# Parameters and initial conditions for the simulation 
p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0] # r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ 
u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ,Nₐ, Sₐ

# Calculate the numerical ODE solution for a 50 days time range
tspan_1 = (0.0, 50.0) 
prob_1 = ODEProblem(single_site!, u0_1, tspan_1, p_1)
sol_1 = solve(prob_1, Tsit5())


"""Abundances plot in time for 50 days"""
plot(sol_1,vars=(0, 1), label="Nⱼ")
plot!(sol_1,vars=(0, 2), label="Nₐ")
title!("1st definition: 50d")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")



# Calculate the numerical ODE solution for a 365 days (one years) time range
tspan_2 = (0.0, 365.0*1) 
prob_2 = ODEProblem(single_site!, u0_1, tspan_2, p_1)
sol_2 = solve(prob_2, Tsit5())


"""Abundances plot in time for 365 days (One Year)"""
plot(sol_2,vars=(0, 1), label="Nⱼ")
plot!(sol_2,vars=(0, 2), label="Nₐ")
title!("1st definition: 1y")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")


# Calculate the numerical ODE solution for a 730 days (two years) time range
tspan_3 = (0.0, 365.0*2) 
prob_3 = ODEProblem(single_site!, u0_1, tspan_3, p_1)
sol_3 = solve(prob_3, Tsit5())


"""Abundances plot in time for 370 days"""

plot(sol_3,vars=(0, 1), label="Nⱼ", color="blue", ylabel="N (Nº individuals)",xlabel="t (days)")
plot!(sol_3,vars=(0, 2), label="Nₐ", color="red",  ylabel="N (Nº individuals)",xlabel="t (days)")
title!("SLC One Site 1")

plot(sol_3,vars=(0, 3), label="Sₐ", color="green", legend=:right, ylabel="S (mm)", xlabel="t (days)")
title!("SLC Once Site 1")



#Fist def final plots: Exploitation rate and Two years simulation:

# For 2 years
plot(sol_3,vars=(0, 1), label="Nⱼ")
plot!(sol_3,vars=(0, 2), label="Nₐ")
title!("1st definition: 2y")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

# Adult size plot by time
# For 2 years
plot(sol_3,vars=(0,3), label="Sₐ") 
title!("1st definition: 2y")
xlabel!("t (days)")
ylabel!("S (mm)")

#=
 - Exploitation of 50% of the adult population:
 During the fist days in the simulation, populations seems to grow and stabilized (fist 50d)
 but with time they drecrease continuously until extintion. 
 
 Posible issue : over-exploitation?
 Posible solution: reduce more the proportion of exploited adult population. 
=#

function Et1(t)
  if modf(t)[1] < 0.5  
    return 0.0
  else
    return 0.37        # Proportion of the adult population that is exploited: 37%
    
  end
end

# Parameters and initial conditions for the diferent proportion of adult exploited: E= 0.37 
p_2 = [0.6, 0.06, 0.05, 0.08, Et1, 1e4, 0.2, 40.0] # r, g, dⱼ, dₐ, E=0.37, K, size_growth_rate, sizeₘₐₓ 
u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ,Nₐ, Sₐ

# Calculate the numerical ODE solution for 2 years (730 days) time range 
prob_3_1 = ODEProblem(single_site!, u0_1, tspan_3, p_2)
sol_3_1 = solve(prob_3_1, Tsit5())



"""Contrast plot of different adult populations exploited (E=0.5 vs E= 0.37)"""

plot(sol_3,vars=(0, 1), label="Nⱼ:E=0.5", legend=:right,style=:dashed, colour="red")
plot!(sol_3,vars=(0, 2), label="Nₐ:E=0.5",style=:solid, colour="red")
plot!(sol_3_1,vars=(0, 1), label="Nⱼ:E=0.37",style=:dashed,colour="blue")
plot!(sol_3_1,vars=(0, 2), label="Nₐ:E=0.37",style=:solid, colour="blue")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")


savefig("C://Code/julia/ODE/Numerical_approximation/Figures/Exploitation_diferencies_N_1.png")
#=
 - Exploitation of 37% of the adult population:
 Now the populations of eggs and adults stabilized with time, number of eggs increase until 7000 individuals aproximately, and
 adults individuals increase until 1500 individuals aproximately.
=#

# Adult size plot by time
# For 50 days

plot(sol_3,vars=(0,3), label="Sₐ:E=0.5",colour="red") 
plot!(sol_3_1,vars=(0,3), label="Sₐ:E=0.37",colour="blue") 
title!("Exploitation diferencies")
xlabel!("t (days)")
ylabel!("S (mm)")

savefig("Exploitation_diferencies_S_1.png")

# With different exploited adult proportions de change in an inversely proportional way.





"""
Second definition: 
Abundances ODE and Size ODE are linked, but no reproductive cicle consideration:
- With function SLC_single_site_S!
"""

#Parameters and initial conditions for this simulations will 
p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0] # r, g, dⱼ, dₐ, E=0.50, K, size_growth_rate, sizeₘₐₓ 
p_2 = [0.6, 0.06, 0.05, 0.08, Et1, 1e4, 0.2, 40.0] # r, g, dⱼ, dₐ, E=0.37, K, size_growth_rate, sizeₘₐₓ 

tspan_3 = (0.0, 365.0*2) ## Two years

# Calculate the numerical ODE solution for a 730 days E:0.5
prob_4 = ODEProblem(single_site_S!, u0_1, tspan_3, p_1) 
sol_4 = solve(prob_4, Tsit5())

# Calculate the numerical ODE solution for a 730 days E:0.37
prob_5 = ODEProblem(single_site_S!, u0_1, tspan_3, p_2) 
sol_5 = solve(prob_5, Tsit5())

#Abundances plot in time 

#For 2 year
plot(sol_4,vars=(0, 1), label = "Nⱼ E:0.5", palette = :darktest)
plot!(sol_4,vars=(0, 2), label= "Nₐ E:0.5", palette = :darktest)
plot!(sol_5,vars=(0, 1), label= "Nⱼ E:0.37", palette = :darktest)
plot!(sol_5,vars=(0, 2), label= "Nₐ E:0.37", palette = :darktest)
title!("Exploitation differencies 2")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")
savefig("Exploitation_diferencies_N_2.png")
#= 
 The abundances presents more fluctuations when the size is linked with eggs abundance, 
 but after 150 days abundances start to decrease until extintion by E:0.37 and no with E:0.5... Why❓
=#

#Adult size plot by time

# For 2 years
plot(sol_4, vars=(0,3), label="E: 0.5")
plot!(sol_5, vars=(0,3), label="E: 0.37")
title!("Adult Size")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

savefig("Exploitation_diferencies_S_2.png")

#= 
 Sizes decrease from 40mm until 30mm aproximately and have more variability in the sizes over time.
=#

"""
The third definition: 
Abundances ODE and Size ODE are linked and have the reproductive cicle consideration:
- With function SLC_single_site_S_X!
"""

#=
 Considering the Exp is a function of t and rate. The ranges used for the exploitation ratis
 in the main.jl "Line 204: exploitation_rates = rand(0.001:0.001:0.003, nsites)", so I decided 
 to choose the max rate (0.003) to observe the most agressive exploitation rate on the population
 in the single site.
=#

# Reproductive cycle
function reproductive_cycle(t)
    if (t % 365) / 365 >= 0.42
      return 1.0
    else
      return 0.0
    end
end

#Exploitation 
# The exploitation cycle is "X·(1-E)" in the ecuations, the exploitation max is E:42
 
function exploit(t, rate)
    if (t % 365) / 365 < 0.37
      return rate
    else
      # return 0.0
      return 0.01
    end
end


# Parameters and initial conditions for the third simulation
#      r,   g,    dⱼ,   dₐ,   E,   K,   size_growth_rate, sizeₘₐₓ, X, rate
p_3 = (0.6, 0.06, 0.05, 0.08, exploit, 1e4, 0.2, 40.0, reproductive_cycle, 0.003)
u0_1 = [1e3, 1e3, 40.0]  ## Nⱼ,Nₐ, Sₐ

tpan_2 =(0.0, 365.0*1)
tpan_3 =(0.0, 365.0*2)

# Calculate the numerical ODE solution for a 365 days (1 year) time range
prob_6 = ODEProblem(single_site_S_X!, u0_1, tspan_2, p_3) 
sol_6 = solve(prob_6, Tsit5())

# Calculate the numerical ODE solution for a 730 days (2 year) time range
prob_7 = ODEProblem(single_site_S_X!, u0_1, tspan_3, p_3) 
sol_7 = solve(prob_7, Tsit5())


#Abundances plot in time 
# For 1 year
plot(sol_6,vars=(0, 1), label="Nⱼ:3º")
plot!(sol_6,vars=(0, 2), label="Nₐ:3º")
title!("3rd definition: 1y")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")


#=
 The tendency present the follow behaviour: extintion perior at the beginning of the year 
 (Na goes to 0), reproduction period at the beginning of the year (Na start to increase 
 and stabilized at 10^4 individuals).
=#

# For 2 years
plot!(sol_7,vars=(0, 1), label="Nⱼ:3º")
plot!(sol_7,vars=(0, 2), label="Nₐ:3º")
title!("3rd definition: 2y")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")


#=
 The tendency present a periodical behaviour.
=#

#Adult size plot by time
plot(sol_7,vars=(0,3), label="Sₐ:3º")
xlabel!("t (days)")
ylabel!("S (mm)")


#Plot for three definition to compare.
#Abundances
plot(sol_3,vars=(0, 1), label="Nⱼ SLC-OS1", style=:dashed, color="red")
plot!(sol_3,vars=(0, 2), label="Nₐ SLC-OS1",style=:solid, color="red")
plot!(sol_5,vars=(0, 1), label="Nⱼ SLC-OS2",style=:dashed, color="blue")
plot!(sol_5,vars=(0, 2), label="Nₐ SLC-OS2",style=:solid, color="blue")
plot!(sol_7,vars=(0, 1), label="Nⱼ SLC-OS3",style=:dashed , color="green")
plot!(sol_7,vars=(0, 2), label="Nₐ SLC-OS3",style=:solid, color="green")
xlabel!("t (days)")
ylabel!("N (Nº individuals)")

# savefig("Model_formulation_SLC_OS_N.png")

#Adult Size
plot(sol_3,vars=(0, 3), label="SLC-OS1",  color="red")
plot!(sol_5,vars=(0, 3), label="SLC-OS2", color="blue")
plot!(sol_7,vars=(0, 3), label="SLC-OS3", color="green")
xlabel!("t (days)")
ylabel!("Sₐ (mm)")
savefig("Model_formulation_SLC_OS_S.png")





### ------------------------------------------------------------------------------
### 2. Numerical aproximation of SLC one site and no reproductive cicle inclusion
### ------------------------------------------------------------------------------


""" ------------------------------------------------------------------------------
-                                    1st definition                              -
---------------------------------------------------------------------------------- """

Exp_lim = 0.9999                 # Exploitation max limit 
m=0.0559                         # Interval of exploitation values 
Expl= 0:m:Exp_lim                # Expoitation values for plotting
tspan = (0.0, 100)                # Time value 
u0 = [1e3,1e2,40]                # Initial conditions of N_e, N_a, S_a

N_et_1 = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
N_at_1 = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
S_at_1 = zeros(Float64,size(Expl)) # Void vector to array the size of adults for diferent exploitation values
c = 0                              # C is the position of the vector N_et, N_at and S_at

for n = 0:m:Exp_lim
  function Et(t)
  if modf(t)[1] < 0.5         # 50% of the adult population is exploited
    return 0.0
  else
    return n                    #For an exploitation value equal to 1, the mathematical result is erratic because the size equation will present a denominator division equal to 0
  end
 end
 p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40] # r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ
 prob_1 = ODEProblem(single_site!, u0, tspan, p_1)
 sol_1 = solve(prob_1, Tsit5())
 c=c+1

 N_et_1[c,] = sol_1[1,end]
 N_at_1[c,] = sol_1[2,end]
 S_at_1[c,] = sol_1[3,end]

end

#Adult size plot by exploitation
plot(Expl,N_et_1,label="Nⱼ")
plot!(Expl,N_at_1,label="Nₐ")
xlims!(0.0,1)
xlabel!("E")
ylabel!("N (nº individuals)")

plot(Expl,S_at_1)
xlabel!("E")
ylabel!("Sₐ (nº individuals)")


""" ------------------------------------------------------------------------------
-                                    2nd definition                              -
---------------------------------------------------------------------------------- """

Exp_lim = 0.9999                 # Exploitation max limit 
m=0.0559                         # Interval of exploitation values 
Expl= 0:m:Exp_lim                # Expoitation values for plotting
tspan = (0.0,100)              # Time range two years
u0 = [1e3,1e3,40]                # Initial conditions of N_e, N_a, S_a


N_et_2 = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
N_at_2 = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
S_at_2 = zeros(Float64,size(Expl)) # Void vector to array the size of adults for diferent exploitation values
c = 0                              # C is the position of the vector N_et, N_at and S_at

for n = 0:m:Exp_lim
  function Et1(t)
  if modf(t)[1] < 0.5
    return 0.0
  else
    return n                    #For an exploitation value equal to 1, the mathematical result is erratic because the size equation will present a denominator division equal to 0
  end
 end
 p_1 = [0.6, 0.06, 0.05, 0.08, Et, 1e4, 0.2, 40.0]
 prob_2 = ODEProblem(SLC_single_site_S!, u0, tspan, p_1)
 sol_2 = solve(prob_2, Tsit5())
 c=c+1

 N_et_2[c,] = sol_2[1,end]
 N_at_2[c,] = sol_2[2,end]
 S_at_2[c,] = sol_2[3,end]

end


plot(Expl,N_et_2,label="Nⱼ")
plot!(Expl,N_at_2,label="Nₐ")
xlims!(0.0,1)
xlabel!("E")
ylabel!("N (nº individuals)")

# In this case abundances doesent change with the exploitation. 


""" ------------------------------------------------------------------------------
-                                    3rd definition                              -
---------------------------------------------------------------------------------- """

Exp_lim = 0.9999                 # Exploitation max limit 
m = 0.0559                       # Interval of exploitation values 
Expl = 0:m:Exp_lim               # Expoitation values for plotting
tspan = (0.0,365*2)              # Time range two years
u0 = [1e3,1e3,40]                # Initial conditions of N_e, N_a, S_a


N_et_3 = zeros(Float64,size(Expl)) # Void vector to array number of eggs for diferent exploitation values
N_at_3 = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
S_at_3 = zeros(Float64,size(Expl)) # Void vector to array the size of adults for diferent exploitation values
c = 0                              # C is the position of the vector N_et, N_at and S_at

for n = 0:m:Exp_lim
  
 # The exploitation cycle is "X·E" in the ecuations, the exploitation max is E:42
 function Exp(t, rate)
  if (t % 365) / 365 < n
    return rate
  else
    # return 0.0
    return 0.00
  end
 end

 # The reproductive cycle is "(1-X)·E" in the theorical equations:
 function X_rc(t)
  if (t % 365) / 365 >= n
    return 1.0
  else
    return 0.00
  end
 end


 p_1 = [0.6, 0.06, 0.05, 0.08, Exp, 1e4, 0.2, 40.0, X_rc, n] # r, g, dⱼ, dₐ, E, K, size_growth_rate, sizeₘₐₓ, X, rate
 prob_3 = ODEProblem(SLC_single_site_S_X!, u0, tspan, p_1)
 sol_3 = solve(prob_3, Tsit5())
 c=c+1

 N_et_3[c,] = sol_3[1,end]
 N_at_3[c,] = sol_3[2,end]
 S_at_3[c,] = sol_3[3,end]

end


plot(Expl,N_et_2,label="Nⱼ")
plot!(Expl,N_at_2,label="Nₐ")
xlims!(0.0,1)
xlabel!("E")
ylabel!("N (nº individuals)")

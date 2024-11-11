#Packages
using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Symbolics
import ForwardDiff.jacobian
#using GLMakie



#Pkg.add("CairoMakie") 
#Pkg.add("DataFrames") 
#Pkg.add("CSV")  
#Pkg.add("DiffEqParamEstim") 
#Pkg.add("Optim") 
#Pkg.add("GLMakie")
#= 
 Formulation of the simple life cicle for one site:

  Population dynamic:

    Simple Life cicle with only one class of population (only one class): Na + Sa

    dNa/dt = X * r * Ne * R * (K - Na/K) - da * Na - (1 - X) * H * Na - da * Na
    +
    dSa/dt = gamma * Sa * (1 - (Sa / (Smax - Smax * H * X)))

    Simple Life Cicle formulation that presents two clases of adult abundances:
    inmature and mature adults: Na, Nm + Sa

    dNa/dt = X * r * Ne * R * (K - Na/K) - gAM * Na - da * Na
    dNm/dt = gAM * Na * (K - (Na + Nm)/K) - (1 - X) * H * Nm - da * Nm
    +
    dSa/dt = gamma * Sa * (1 - (Sa / (Smax - Smax * H * X)))

 Formulation for complex life cicles for one site: different numbers of clases.

  Population dynamics:

    Complex Life Cicle with 2 clases of population: Ne, Na + Sa
    
    dNe/dt = X * r * Na * R - gEA * Ne - de * Ne
    dNa/dt = gEA * Ne * (K - Na/K) - (1 - X) * H * Na - da * Na
    +
    dSa/dt = gamma * Sa * (1 - (Sa / (Smax - Smax * H * X)))

    Complex Life Cicle with 5 clases of population: Ne, Nt, Nv, Nj, Na + Sa

    dNe/dt = X * r * Na * R - gEA * Ne - de * Ne
    dNt/dt = gET * Ne - gTV * Nt - dt * dNt
    dNv/dt = gTV * Nt * (K - Nv/K) - gVJ * Nv - dv * Nv
    dNj/dt = gVJ * Nv * (K - Nj/K) - gJA * Nj - dj * Nj
    dNa/dt = gJA * Nj * (K - Na/K) - (1 - X) * H * Na - da * Na
    +
    dSa/dt = gamma * Sa * (1 - ( Sa / (Smax - Smax * H * X)))

 where i = nº of species == ["Patella ordinaria", "Patella aspera"]

Parameters and variables:
 - Ne = eggs abundance
 - Nt = trocophore abuncance
 - Nv = veliger abuncance
 - Nj = juvenile abuncance
 - Na = adults abuncance (when Nm is in the equations, Na means non matured adults)
 - Nm = matured adult abundaces
 - Sa = adults size (Average sizes Before MPS+FULL = [46.44,44.45])
 - r = population growth rate [9.17,5.03]
 - R = reproductive capacity
 - K = carrying capacity (k = 1e^4)
 - X = Reproductive period [1,0] 
 - (1-X) = Exploitation periosd [1,0]
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
=#

# Full access scenarios

function SLC!(du, u, p, t)
   Na, Sa = u
   r, K, H, da, Smax, gamma = p

   Saverage = du[2]
   Smaturity = calculate_size_at_first_maturity(Saverage)
   
  du[1] = dNa = X(t) * r * Na * reproduction_capacity(Saverage, Smaturity, Smax) * ((K - Na)/K) - (1 - X(t)) * H * Na - (da * Na) 
  du[2] = dSa = gamma * Sa * (1 - Sa / (Smax - (1 - H)))
end

#=
function aSLC!(du, u, p, t)
  Na,Nm, Sa = u
  r, K, H, X, da, Smax = p
  Saverage = du[3]
  Smaturity = calculate_size_at_first_maturity(Saverage)

 du[1] = dNa = X(t) * r * Na * Rep_cap(Saverage, Smaturity, Smax) * (K - Na/K) - (da * Na) 
 du[2] = dNm = gAM * Na * (K - (Na + Nm)/K) - (1 - X(t)) * H * Nm - (da * Nm) 
 du[3] = dSa = gamma * Sa * (1 - Sa / (Smax - (1 - H)))
end
=#
function SCLC!(du, u, p, t)
  Ne, Na, Sa = u
   r, K, H, gEA, de, da, Smax, gamma = p
   Saverage = du[3]
   Smaturity = calculate_size_at_first_maturity(Saverage)
   
  du[1] = dNe = X(t) * r * Na * reproduction_capacity(Saverage, Smaturity, Smax) - (de * Ne) - (gEA * Ne)
  du[2] = dNa = gEA * Ne * ((K - Na)/K) - (da * Na) - ((1 - X(t))* H * Na)
  du[3] = dSa = gamma * Sa * (1 - Sa / (Smax - (1 * (1 - X(t)) * H)))
end
#=
function aSCLC!(du, u, p, t)
  Ne, Na, Sa = u
   r, K, H, X, gEA,gAM, de, da, Smax, gamma = p
   Saverage = du[4]
   Smaturity = calculate_size_at_first_maturity(Saverage)

  du[1] = dNe = X(t) * r * Na * Rep_cap(Saverage, Smaturity, Smax) - (de * Ne) - (gEA * Ne)
  du[2] = dNa = gEA * Ne * (K - Ne/K) - gAM * Na - da * Na 
  du[3] = dNm = gAM * Na * (K - (Na + Nm)/K) - (1 - X(t)) * H * Nm - (da * Nm) 
  du[4] = dSa = gamma * Sa * (1 - Sa / (Smax - (1 * (1 - X(t)) * H)))
end
=#
function CLC!(du, u, p, t)
  Ne, Nt, Nv, Nj, Na, Sa = u
   r, K, H, g ,de, dt, dv, dj, da, Smax, gamma = p
   Saverage = du[6]
   Smaturity = calculate_size_at_first_maturity(Saverage)

  du[1] = dNe = X(t) * r * Na * reproduction_capacity(Saverage, Smaturity, Smax)- de * Ne - g[1] * Ne
  du[2] = dNt = g[1] * Ne - g[2] * Nt - dt * Nt
  du[3] = dNv = g[2] * Nt * ((K - Nt)/K) - g[3] * Nv - dv * Nv
  du[4] = dNj = g[3] * Nv * ((K - Nj)/K) - g[4] * Nj - dj * Na
  du[5] = dNa = g[4] * Nj * ((K - Na)/ K) - da * Na - (1 - X(t)) * H * Na
  du[6] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - X(t))))
end

#=
function aCLC!(du, u, p, t)
  Ne, Nt, Nv, Nj, Na, Sa = u
   r, K, H, X, gET, gTV, gVJ, gJA, de, dt, dv,dj, da, Smax, gamma = p
   Saverage = du[7]
   Smaturity = calculate_size_at_first_maturity(Saverage)

  du[1] = dNe = X(t) * r * Na * Rep_cap(Saverage, Smaturity, Smax) - de * Ne - gEA * Ne
  du[2] = dNt = gET * Ne - gTV*Nt - dt * Nt
  du[3] = dNv = gTV * Nt * (K - Nt/K) - gVJ * Nv - dv * Nv
  du[4] = dNj = gVJ * Nv * (K - Nj/K) - gJA * Nj - dj * Na
  du[5] = dNa = gJA * Nv * (K - Na/K) - gAM * Na - da * Na 
  du[6] = dNm = gAM * Na * (K - (Na + Nm)/K) - (1 - X(t)) * H * Nm - (da * Nm) 
  du[7] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - X(t))))
end
=#
# Reproductive Cycle (X=1)

function X(t)
  if (t % 365) / 365 >= 0.42
    return 1.0 # Reproductive Cycle
  else
    return 0.0 # Exploitation Cycle
  end
end

function calculate_size_at_first_maturity(current_avg_size)
  M = 1.34 * (current_avg_size) - 28.06
end

"Return the reproduction capacity (between 0 and 1) given the current average size and size at first maturity and maximum size"
reproduction_capacity(Saverage, Smaturity, Smax) = min(max(0.5 * (1.0 + (Saverage - Smaturity) / (Smax - Smaturity)), 0.0), 1.0)


# Population Growth rate estimation (r=reggs):

oocytes_po = 385613.0                # Average: Patella ordinaria (nº of Eggs)
oocytes_pa = 73029.0                   # Average: Patella aspera (nº of Eggs)
oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.

re = reggs / 500     # Population growth rate
gs = [0.998611, 0.971057, 0.4820525, 0.00629]

Kt = 640000.0        # Carrying capacity
rates = [0.639,0.57] # Exploitation rate (H)

gEA = 0.006         # Instant conversion between stages.

# Natural mortality rates:
# see estimate_mortality_rates.jl for how these values were estimated.

de_ = 0.99 / 365
dt_ = 0.717 / 365
dv_ = 0.392 / 365
dj_ = 0.315 / 365 
da_ = 0.1175 / 365
d_= [0.55,0.59]    # Natural mortality rate for speciesd

Sm = 56.0             # Maximum size for adults

gammas = [0.32,0.36] # Adult growth rate

i = [1,2]           # Species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)


p_SLC_po = [re[1], Kt, rates[1], d_[1], Sm, gammas[1]]
p_SLC_pa = [re[2], Kt, rates[2], d_[2], Sm, gammas[2]]


p_SCLC_po = [re[1], Kt, rates[1], gEA, de_, da_, Sm, gammas[1]]
p_SCLC_pa = [re[2], Kt, rates[2], gEA, de_, da_, Sm, gammas[2]]



p_CLC_po = [re[1], Kt, rates[1], gs, de_, dt_, dv_, dj_, da_, Sm, gammas[1]]
p_CLC_pa = [re[2], Kt, rates[2], gs, de_, dt_, dv_, dj_, da_, Sm, gammas[2]]

start_year = 2006
end_year = 2018

nyears = end_year - start_year + 1
ndays = nyears * 365.0
t_span = (0.0, ndays) # Temporal ranges for simulations.
 
u0_SLC_po_full = [1e4, 43.41]    # Patella ordinaria 
u0_SLC_pa_full = [1e4, 45.72]    # Patella aspera

u0_SCLC_po_full = [1e4, 1e4, 43.41]    # Patella ordinaria 
u0_SCLC_pa_full = [1e4, 1e4, 45.72]    # Patella aspera


u0_CLC_po_full = [1e4, 1e4, 1e4, 1e4, 1e4, 43.41]    # Patella ordinaria 
u0_CLC_pa_full = [1e4, 1e4, 1e4, 1e4, 1e4, 45.72]    # Patella aspera



#SLC
prob_SLC_full = ODEProblem(SLC!, u0_SLC_po_full, t_span, p_SLC_po) 
sol_SLC_full = solve(prob_SLC_full, Tsit5())


fig1 = Figure()
ax1 = Axis(fig1[1, 1])
lines!(ax1, sol_SLC_full.t, [u[1] for u in sol_SLC_full.u], yscale=:log10, label="Na")
#save("SLC_N_Full_access.png", fig1, dpi = 300)

fig11 = Figure()
ax11 = Axis(fig11[1, 1])
lines!(ax11, sol_SLC_full.t, [u[2] for u in sol_SLC_full.u], yscale=:log10, label="Sa")
save("SLC_S_Full_access.png", fig11, dpi = 300)






# SCLC
prob_SCLC_full = ODEProblem(SCLC!, u0_SCLC_po_full, t_span, p_SCLC_po) 
sol_SCLC_full = solve(prob_SCLC_full, Tsit5())

fig2 = Figure()
ax2 = Axis(fig2[1, 1])
lines!(ax2, sol_SCLC_full.t, [u[1] for u in sol_SCLC_full.u], yscale=:log10, label="Ne")
lines!(ax2, sol_SCLC_full.t, [u[2] for u in sol_SCLC_full.u], yscale=:log10, label="Na")
save("SCLC_N_Full_access.png", fig2, dpi = 300)

fig3 = Figure()
ax3 = Axis(fig3[1, 1])
lines!(ax3, sol_SCLC_full.t, [u[3] for u in sol_SCLC_full.u], yscale=:log10, label="Sa",)
save("SCLC_S_Full_access.png", fig3, dpi = 300)




#CLC
prob_CLC_full = ODEProblem(CLC!, u0_CLC_po_full, t_span, p_CLC_po) 
sol_CLC_full = solve(prob_CLC_full, Tsit5())



fig4 = Figure()
ax4 = Axis(fig4[1, 1])
lines!(ax4, sol_CLC_full.t, [u[1] for u in sol_CLC_full.u], yscale=:log10, label="Ne",
 xlabel = "time", ylabel = "Na (nº of individuals)")
lines!(ax4, sol_CLC_full.t, [u[2] for u in sol_CLC_full.u], yscale=:log10, label="Nt")
lines!(ax4, sol_CLC_full.t, [u[3] for u in sol_CLC_full.u], yscale=:log10, label="Nv")
lines!(ax4, sol_CLC_full.t, [u[4] for u in sol_CLC_full.u], yscale=:log10, label="Nj")
lines!(ax4, sol_CLC_full.t, [u[5] for u in sol_CLC_full.u], yscale=:log10, label="Na")
save("CLC_N_Full_access.png", fig4, dpi = 300)


fi5 = Figure()
ax5 = Axis(fig5[1, 1])
lines!(ax5, sol_CLC_full.t, [u[6] for u in sol_CLC_full.u], yscale=:log10, label="Sa",
xlabel = "time", ylabel = "Sa (mm)")
save("CLC_S_Full_access.png", fig4, dpi = 300)

#=
Analytical Aproximation of SLC in One Place.

function SLC!(du, u, p, t)
  Na, Sa = u
  i, r, K, H, X, da, Smax, gamma = p

  Saverage = du[2]
  Smaturity = calculate_size_at_first_maturity(Saverage)

 du[1] = dNa = X(t) * r * Na * Rep_cap(Saverage, Smaturity, Smax) * (K - Na/K) - (1 - X(t)) * H(i) * Na - (da * Na) 
 du[2] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 - H[i])))
end
=#

Exp_lim = 0.9999                 # Exploitation max limit 
m = 0.0001                       # Interval of exploitation values 
Expl = 0:m:Exp_lim               # Expoitation values for plotting
tspan = (0.0,365*2)              # Time range two years
K = 640000.0    # Capacidad de carga                
# Initial conditions of N_a, S_a

N_at = zeros(Float64,size(Expl)) # Void vector to array number of adults for diferent exploitation values
S_at = zeros(Float64,size(Expl)) # Void vector to array the size of adults for diferent exploitation values
c = 0                              # C is the position of the vector N_et, N_at and S_at

for H in 0:m:Exp_lim
  
 Smax = 56
 X_val = 1
 Saverage = mean([42 56])
 Smaturity = 1.34 * Saverage - 28.06
 R = min(max(0.5 * (1.0 + (Saverage - Smaturity) / (Smax - Smaturity)), 0.0), 1.0)

 Na_values = (K/2)*(K+(-d+H*(1 - X_val))/(X_val*r*R))
 
 Sa_values = (Smax*(1-H))/2

 
 c=c+1
 N_at[c,] = Na_values
 S_at[c,] = Sa_values
 
end


lines(Expl,N_at,label="Nₐ")
xlims!(0.0,1)
xlabel!("E")
ylabel!("N (nº individuals)")
savefig!("SLC_N_prima.png")


lines(Expl,S_at,label="Sₐ")
xlims!(0.0,1)
xlabel!("E")
ylabel!("S(mm)")
savefig!("SLC_S_prima.png")



#Check SLCmeta.jl

function my_ode!(du, u, t, p)
    x1, x2, y1, y2 = u
    re1, re2, K, H1, H2, X, R1, R2, g1, g2, d1, d2, c12, c21, Smax = p
    
    #Parameter values species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)
    # Population Growth rate estimation (r=reggs):

   # oocytes_po = 385613                  # Average: Patella ordinaria (nº of Eggs)
   # oocytes_pa = 73029                   # Average: Patella aspera (nº of Eggs)
   # oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
   # reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.
   # re = reggs / 500     # Population growth rate: re[1,1] re[2,1]
    
    du[1] = dx1 = (X(t) * re1 * R1 * x1)*((K - x1)/K) - (d1 * x1) - H1*(1 - X(t))*x1 - c12*x1*x2
    du[2] = dx2 = (X(t) * re2 * R2 * x2)*((K - x2)/K) - (d2 * x2) - H2*(1 - X(t))*x2 - c21*x1*x2
    du[3] = dy1 = (g1 * y1) * (1 - (y1/(Smax * (1 - H1 * (1-X(t))))))
    du[4] = dy2 = (g2 * y2) * (1 - (y2/(Smax * (1 - H2 * (1-X(t))))))
end

# Reprouctive cycle X
function rep(t)
  if (t % 365) / 365 >= 0.42
    return 1.0 # Reproductive Cycle
  else
    return 0.0 # Exploitation Cycle
  end
end

#Initial parameters
u0 = [100.0, 100.0, 25.0, 25.0]
tspan = (0.0, 1000.0)

R1 = 1
R2 = 1
K = 640000          # Carrying capacity 
H1 = 0.639
H2 = 0.57 # Exploitation rate (H)
g1 = 0.32
g2 = 0.36 # Adult growth rate
 #rates2 = [0.02,0.01]
 #gEA = 0.006          # Instant conversion between stages.
d1 = 0.55
d2 = 0.59  # Natural mortality rate for adults
Smax = 56              # Maximum size for adults
re1 = 0.32
re2 = 0.36 # Adult growth rate
c12 = 0.05 #competition term species 2 on 1
c21 = 0.05 #competition term species 1 on 2

# Exploitation gradient
n_values = 0:0.05:1

# Cubic matrix for save simulation outputs 
du_t_n = zeros(length(u0), length(tspan[1]:tspan[2]), length(n_values))
sol = zeros(length(u0), length(tspan[1]:tspan[2]))
# Simulation loop for each exploitation value of the gradient

#for (i, n) in enumerate(n_values)
 n_values = 0.5

    p_solve = re1, re2, K, n_values, H2, rep, R1, R2, g1, g2, d1, d2, c12, c21, Smax
    # Define ODE Problem with a new H1 value
    prob = ODEProblem(my_ode!, u0, tspan, p_solve)
    
    # Solve ODE problem and save solution in a cubic matrix
    sol = solve(prob, Tsit5()) # There is an error when I run this line and I don't know how to fix it.
    du_t_n[:, :, i] = sol
#end

fig_acum = plot!(du_t_n[1,500,:], xlabel="Exploitation", ylabel="State Variables", label=["n"])


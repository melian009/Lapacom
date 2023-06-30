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
using Plots


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
   i, r, K, H, X, da, Smax, gamma = p

   Saverage = du[2]
   Smaturity = calculate_size_at_first_maturity(Saverage)

  du[1] = dNa = X(t) * r * Na * Rep_cap(Saverage, Smaturity, Smax) * ((K - Na)/K) - (1 - X(t)) * H(i) * Na - (da * Na) 
  du[2] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 - H[i])))
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
   i,r, K, H, X, gEA, de, da, Smax, gamma = p
   Saverage = du[4]
   Smaturity = calculate_size_at_first_maturity(Saverage)

  du[1] = dNe = X(t) * r[i] * Na * Rep_cap(Saverage, Smaturity, Smax) - (de * Ne) - (gEA * Ne)
  du[2] = dNa = gEA * Ne * ((K - Na)/K) - (da * Na) - ((1 - X(t))* H[i] * Na)
  du[3] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 * (1 - X(t)) * H[i])))
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
  i, r, K, H, X, g ,de, dt, dv, dj, da, Smax, gamma = p
   Saverage = du[6]
   Smaturity = calculate_size_at_first_maturity(Saverage)

  du[1] = dNe = X(t) * r[i] * Na * Rep_cap(Saverage,Smaturity,Smax) - de * Ne - g[1] * Ne
  du[2] = dNt = g[1] * Ne - g[2] * Nt - dt * Nt
  du[3] = dNv = g[2] * Nt * ((K - Nt)/K) - g[3] * Nv - dv * Nv
  du[4] = dNj = g[3] * Nv * ((K - Nj)/K) - g[4] * Nj - dj * Na
  du[5] = dNa = g[4] * Nj * ((K - Na )/ K) - da * Na - (1 - X(t)) * H[i] * Na
  du[6] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - Smax * H * (1 - X(t))))
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

function rep(t)
  if (t % 365) / 365 >= 0.42
    return 1.0 # Reproductive Cycle
  else
    return 0.0 # Exploitation Cycle
  end
end

function calculate_size_at_first_maturity(Sav)
  M = 1.34 * (Sav) - 28.06
end

# "Return the reproduction capacity (between 0 and 1) given the current average size
#  and size at first maturity and maximum size".
Rep_cap(Saverage, Smaturity, Smax) = min(max(0.5 * (1.0 + (Saverage - Smaturity) / (Smax - Smaturity)), 0.0), 1.0)


# Population Growth rate estimation (r=reggs):

oocytes_po = 385613                  # Average: Patella ordinaria (nº of Eggs)
oocytes_pa = 73029                   # Average: Patella aspera (nº of Eggs)
oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.

re = reggs / 500     # Population growth rate
gs = [0.998611, 0.971057, 0.4820525, 0.00629]

Kt = 640000          # Carrying capacity
rates = [0.639,0.57] # Exploitation rate (H)

gEA = 0.006          # Instant conversion between stages.

# Natural mortality rates:
# see estimate_mortality_rates.jl for how these values were estimated.

de_ = 0.99 / 365
dt_ = 0.717 / 365
dv_ = 0.392 / 365
dj_ = 0.315 / 365 
da_ = 0.1175 / 365
d_= [0.55,0.59]    # Natural mortality rate for speciesd

Sm = 56              # Maximum size for adults

gammas = [0.32,0.36] # Adult growth rate

i = [1,2]            # Species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)


p_SLC_po = [i[1], re[1], Kt, rates[1], rep, da_, Sm, gammas[1]]
p_SLC_pa = [i[2], re[2], Kt, rates[2], rep, da_, Sm, gammas[2]]


p_SCLC_po = [i[1],re[1], Kt, rates[1], rep, gEA, de_, da_, Sm, gammas[1]]
p_SCLC_pa = [i[2],re[2], Kt, rates[2], rep, gEA, de_, da_, Sm, gammas[2]]


p_CLC_po = [i[1],re[1], Kt, rates[1], rep, gs, de_, dt_, dv_, dj_, da_, Sm, gammas[1]]
p_CLC_pa = [i[2],re[1], Kt, rates[2], rep, gs, de_, dt_, dv_, dj_, da_, Sm, gammas[2]]

t_span= (1,30000) # Temporal ranges for simulations: 2 years.
 
u0_SLC_po_full = [1e4, 43.41]    # Patella ordinaria 
u0_SLC_pa_full = [1e4, 45.72]    # Patella aspera

u0_SCLC_po_full = [1e4, 1e4, 43.41]    # Patella ordinaria 
u0_SCLC_pa_full = [1e4, 1e4, 45.72]    # Patella aspera


u0_CLC_po_full = [1e4, 1e4, 1e4, 1e4, 1e4, 43.41]    # Patella ordinaria 
u0_CLC_pa_full = [1e4, 1e4, 1e4, 1e4, 1e4, 45.72]    # Patella aspera


prob_SLC_full = ODEProblem(SLC!, u0_SLC_po_full, t_span, p_SLC_po) 
sol_SLC_full = solve(prob_SLC_full, Tsit5())

prob_SCLC_full = ODEProblem(SCLC!, u0_SCLC_po_full, t_span, p_SCLC_po) 
sol_SCLC_full = solve(prob_SCLC_full, Tsit5())

prob_CLC_full = ODEProblem(CLC!, u0_CLC_po_full, t_span, p_CLC_po) 
sol_CLC_full = solve(prob_CLC_full, Tsit5())


# prob_pa_full = ODEProblem(CLC!, u0_pa_full, t_span, p_po) 
# sol_pa_full = solve(prob_pa_full, Tsit5())


Plots.plot(sol_CLC_full, vars=(0,1),  label= "Ne (Full access)")
Plots.plot!(sol_CLC_full, vars=(0,2), label= "Nt (Full access)")
Plots.plot!(sol_CLC_full, vars=(0,3), label= "Nv (Full access)")
Plots.plot!(sol_CLC_full, vars=(0,4), label= "Nj (Full access)")
Plots.plot!(sol_CLC_full, vars=(0,5), label= "Na (Full access)")
Plots.title!("'Patella ordinaria'")
Plots.xlabel!("t (days)")
Plots.ylabel!("LOG10(N) (Nº individuals)")
#savefig!("CLC_SS_po_N_Full_access_log.png")


plot(sol_CLC_full, vars=(0,6),  label= "Ne (Full access)")
title!("'Patella ordinaria'")
xlabel!("t (days)")
Plots.ylabel!("N (Nº individuals)")
savefig!("CLC_SS_po_N_Full_access.png")


#=
Plots.plot(sol_pa_full, vars=(0,1), yscale=:log10,  label= "Ne (Full access)")
Plots.plot!(sol_pa_full, vars=(0,2), yscale=:log10, label= "Nt (Full access)")
Plots.plot!(sol_pa_full, vars=(0,3), yscale=:log10, label= "Nv (Full access)")
Plots.plot!(sol_pa_full, vars=(0,4), yscale=:log10, label= "Nj (Full access)")
Plots.plot!(sol_pa_full, vars=(0,5), yscale=:log10, label= "Na (Full access)")
Plots.title!("'aatella ordinaria'")
Plots.xlabel!("t (days)")
Plots.ylabel!("LOG10(N) (Nº individuals)")
savefig!("CLC_SS_pa_N_Full_access.png")
=#

```
Analytical Aproximation of SLC in One Place.

function SLC!(du, u, p, t)
  Na, Sa = u
  i, r, K, H, X, da, Smax, gamma = p

  Saverage = du[2]
  Smaturity = calculate_size_at_first_maturity(Saverage)

 du[1] = dNa = X(t) * r * Na * Rep_cap(Saverage, Smaturity, Smax) * (K - Na/K) - (1 - X(t)) * H(i) * Na - (da * Na) 
 du[2] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 - H[i])))
end
```

function X(t)
  if (t % 365) / 365 >= 0.42
      return 1.0 # Reproductive Cycle
  else
      return 0.0 # Exploitation Cycle
  end
end

function simulate_NS_values()
  No = 10000.0    # Initial abuncance
  So = 43.41      # Initial size
  r = 2515.4/500    # Reproductve rate
  K = 640000.0    # Carrying capacity
  d = 0.55        # Natual mortality
  Smax = 56.0     # Mamimum size
  gamma = 0.34    # Growth rate

  t_max = 365*2    # End time for the simulation
  step_size = 1   # Time intervales (days)

  H_values = 0.0:0.05:1.0  # Exploitatation rates

  NS_matrix = zeros(length(H_values), 3)   # Output matrix for N, S and H estimated values

  for (i, H) in enumerate(H_values)
      t_values = 0:step_size:t_max
      Na_values = zeros(length(t_values))
      Sa_values = zeros(length(t_values))

      for (j, t) in enumerate(t_values)
          X_val = X(t)
          Saverage = mean(Sa_values[1:j])
          Smaturity = 1.34 * Saverage - 28.06
          R = min(max(0.5 * (1.0 + (Saverage - Smaturity) / (Smax - Smaturity)), 0.0), 1.0)
          #Stability point (N=0,S=0)
          #Na_values[j] = No * exp((X_val * r * R * K - d - H * (1 - X_val) * t))
          #Sa_values[j] = So * exp(gamma * t)
          #Stability Point (N=(X*R*r*K^2-d-H*(1-x)),S=(Smax*(1-H*(1-X))))
          Na_values[j] = No * exp(X_val * r * R * K * (1 - 2 * X_val * r * R) + ((2 * X_val * r * R)/K - 1)*(d + H * (1 - X_val))) + X_val * r * R * (K^2) - d - H * (1 - X)
          Sa_values[j] = Smax * (1 - H * (1-X_val)) - So * exp(gamma * t)
      end

      N_median = median(Na_values)
      S_median = median(Sa_values)

      NS_matrix[i, 1] = H
      NS_matrix[i, 2] = N_median
      NS_matrix[i, 3] = S_median

  end

  return NS_matrix
end

NS_matrix = simulate_NS_values()

# Plot S vs H
plot(NS_matrix[:, 1], NS_matrix[:, 3], xlabel = "H", ylabel = "S", label = "S vs H", legend=:topleft)
title!("S values for different H")

# Plot N vs H
plot(NS_matrix[:, 1], NS_matrix[:, 2], xlabel = "H", ylabel = "N", label = "N vs H", legend=:topleft)
title!("N values for different H")
#Packages
using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using GlobalSensitivity
#using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Symbolics
import ForwardDiff.jacobian
import Plots

function SLC!(du, u, p, t)
   Na, Sa = u
   i, r, K, H, X, da, Smax, gamma = p
   Saverage = du[2]
   Smaturity = calculate_size_at_first_maturity(Saverage)

  du[1] = dNa = X(t) * r * Na * Rep_cap(Saverage, Smaturity, Smax) * (K - Na/K) - (1 - X(t)) * H(i) * Na - (da[i] * Na) 
  du[2] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 - H[i])))
end


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


t_span= (1,3000) # Temporal ranges for simulations: 2 years.
 
u0_SLC_po_full = [1e4, 43.41]    # Patella ordinaria 
u0_SLC_pa_full = [1e4, 45.72]    # Patella aspera

prob_SLC_full = ODEProblem(SLC!, u0_SLC_po_full, t_span, p_SLC_po) 
sol_SLC_full = solve(prob_SLC_full, Tsit5())


#Plots.plot(sol_po_full, vars=(0,1), yscale=:log10,  label= "Ne (Full access)")
#Plots.plot!(sol_po_full, vars=(0,2), yscale=:log10, label= "Nt (Full access)")
#Plots.title!("'Patella ordinaria'")
#Plots.xlabel!("t (days)")
#Plots.ylabel!("LOG10(N) (Nº individuals)")
#savefig!("CLC_SS_po_N_Full_access_log.png")



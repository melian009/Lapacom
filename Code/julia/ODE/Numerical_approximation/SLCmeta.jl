using Pkg
Pkg.add("DifferentialEquations")
function my_ode!(du, u, t)
    i, r, K, H, X, g, da, Smax, gamma = p
    du[1] = (X * r1 * R1 * x1)*((K - x1)/K) - (d1 * x1) - H*(1 - X)*x1 - c12*x1*x2
    du[2] = (X * r2 * R2 * x2)*((K - x2)/K) - (d2 * x2) - H*(1 - X)*x2 - c21*x1*x2
    du[3] = (g1 * y1) * (1 - (y1/(Smax1 * (1 - H * (1-X)))))
    du[4] = (g2 * y2) * (1 - (y2/(Smax2 * (1 - H * (1-X)))))
end


unction SLC_metapop_before!(du, u, p, t)
  Na, Sa = u
   i, r, K, H, X, g, da, Smax, gamma = p
  #du[1] = dNe = (X(t) * r[i] * Na * ((K - Na) / K)) - (de[i] * Ne) - (g * Ne)
  du[2] = dNa = (g * Ne * ((K - Ne) / K)) - (da[i] * Na) - (H[i] * Na)
  du[3] = dSa = gamma[i] * Sa * (1 - Sa / (Smax - (1 - H[i])))
end

u0 = [100.0, 100.0, 25, 25]  # Initial conditions
tspan = (0.0, 1000.0)  # Time span for the simulation (from t=0 to t=1000)

using DifferentialEquations
solver = Tsit5()
prob = ODEProblem(my_ode!, u0, tspan)
sol = solve(prob, solver)

using Plots
plot(sol, xlabel="Time", ylabel="State Variables", label=["x1" "x2" "y1" "y2"])


# Population Growth rate estimation (r=reggs):

oocytes_po = 385613                  # Average: Patella ordinaria (nº of Eggs)
oocytes_pa = 73029                   # Average: Patella aspera (nº of Eggs)
oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
reggs = oocytes / (365 * 0.42)       # conversion rate of adults to eggs.

re = reggs / 500     # Population growth rate
K = 640000          # Carrying capacity
rates = [0.639,0.57] # Exploitation rate (H)
rates2 = [0.02,0.01]
gEA = 0.006          # Instant conversion between stages.
da_ = [0.55,0.59]    # Natural mortality rate for adults
Sm = 56              # Maximum size for adults
gammas = [0.32,0.36] # Adult growth rate
i = [1,2]            # Species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)


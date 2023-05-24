using Optim
using DifferentialEquations

S0 = 48 # initial size
S1 = S0 + 2.0
t = 730.0

function sizeonly_gompertz!(du, u, p, t)
  E = 0.0 # You might want to adjust this based on your specific needs
  sizeₘₐₓ = S1  # maximum size
  du[1] = p * u[1] * exp(-u[1] / (sizeₘₐₓ - (sizeₘₐₓ * E)))
end

function size_growth_rate(r)
  u0 = [S0]
  tspan = (0.0, t)
  prob = ODEProblem(sizeonly_gompertz!, u0, tspan, r)
  sol = solve(prob)
  (S1 - sol[end][1])^2
end

opt = optimize(size_growth_rate, 0.0, 1.0)

r_opt = opt.minimizer  # 0.00014898749263737575

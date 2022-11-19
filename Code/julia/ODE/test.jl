using Pkg
Pkg.activate(".")
using OrdinaryDiffEq
using CairoMakie


# u0 = [[1_800.0 for j in 1:2] for i in 1:8]
# u0 = reduce(vcat, u0)
u0 = reshape([1_800.0 for j in 1:24], 8, 3)
r = [10.0, 0.99, 0.92]
d = [0.001, 0.001, 0.001]
K = 40_000
α = [0.1, 0.1, 0.1]
p = [r, d, K, α]

function t1!(du, u, p, t)
  counter = 0
  for site in 1:8
    counter += 1
    stage = 1
    prev_stage = 2
    du[counter] = (p[1][stage] * u[site, prev_stage] * ((p[3] - u[site, prev_stage]) / p[3])) -
                  (p[1][stage+1] * u[site, stage]) -
                  (p[2][stage] * u[site, stage])
  end
end


## defining model parameters
function exploit(t, rate)
  if (t % 365) / 365 < 0.42
    return rate
  else
    return 0.0
  end
end

function reproductive_cycle(t)
  if (t % 365) / 365 >= 0.42
    return 1.0
  else
    return 0.0
  end
end

function t2!(du, u, p, t)
  nsites = 8
  nstages = 3

  counter = 0
  for site in 1:nsites

    #stage 1
    counter += 1
    stage = 1
    prev_stage = nstages
    dispersal_probs1 = fill(0.1, nsites)
    dispersal_probs2 = fill(0.09, nsites)
    du[counter] = (reproductive_cycle(t) * p[1][stage] * u[site, prev_stage] * ((p[3] - u[site, prev_stage]) / p[3])) -
                  (p[1][stage+1] * u[site, stage]) -
                  (p[2][stage] * u[site, stage]) +
                  (sum(dispersal_probs1 .* u[:, stage])) -
                  (u[site, stage] * sum(dispersal_probs2))
    
    counter += 1
    stage = 2
    prev_stage = 1
    dispersal_probs1 = fill(0.1, nsites)
    dispersal_probs2 = fill(0.09, nsites)
    du[counter] = (p[1][stage] * u[site, prev_stage]) -
                  (p[1][stage+1] * u[site, stage]) -
                  (p[2][stage] * u[site, stage]) +
                  (sum(dispersal_probs1 .* u[:, stage])) -
                  (u[site, stage] * sum(dispersal_probs2))

    counter += 1
    stage = 3
    prev_stage = 2
    dispersal_probs1 = fill(0.1, nsites)
    dispersal_probs2 = fill(0.09, nsites)
    du[counter] = (p[1][stage] * u[site, prev_stage]) -
              (p[1][stage] * u[site, stage]) -
              (p[2][stage] * u[site, stage]) +
              (sum(dispersal_probs1 .* u[:, stage])) -
              (u[site, stage] * sum(dispersal_probs2))
              
  end
end

tspan = (100.0, 200.0)
prob = ODEProblem(t2!, u0, tspan, p)
sol = solve(prob, Rosenbrock23(), dt=0.0001, adaptive=false)

fig = Figure()
ax1 = Axis(fig[1, 1])
stage = 1
lines!(ax1, sol.t, [sol.u[i][1, stage] for i in 1:length(sol.t)], yscale=:log10)
for site in 2:8
  lines!(ax1, sol.t, [sol.u[i][site, stage] for i in 1:length(sol.t)])
end
save("figs/test_stage=$stage.pdf", fig)
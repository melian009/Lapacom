using Pkg
Pkg.add("DifferentialEquations")
function my_ode!(du, u, p, t)
    x, y = u
    du[1] = x^2 - y
    du[2] = x - y^2
end

u0 = [1.0, 2.0]  # Initial conditions
tspan = (0.0, 10.0)  # Time span for the simulation (from t=0 to t=10)

using DifferentialEquations
solver = Tsit5()
prob = ODEProblem(my_ode!, u0, tspan)
sol = solve(prob, solver)

using Plots
plot(sol, xlabel="Time", ylabel="State Variables", label=["x" "y"])


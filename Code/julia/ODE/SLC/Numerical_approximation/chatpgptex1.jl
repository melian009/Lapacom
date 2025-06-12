using DifferentialEquations
using Plots


function ode_system!(du, u, p, t)
    r = p
    r = -0.5
    du[1] = r * u[1] + u[1]^2 - u[1] * u[2]
    du[2] = u[1] - u[1] * u[2]
end


u0 = [1.0, 1.0]  # Initial conditions for [x, y]
tspan = (0.0, 10.0)  # Time span for the simulation (from t=0 to t=10)



# Define the ODE problem
prob = ODEProblem(ode_system!, u0, tspan)

# Solve the ODE numerically
sol = solve(prob)

# Plot the solution
plot(sol, xlabel="Time", ylabel="State Variables", label=["x" "y"])


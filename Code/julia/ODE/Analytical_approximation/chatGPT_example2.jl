using SymPy

# Define the symbolic variables
t, y = symbols("t y")

# Define the nonlinear ODE, for example: dy/dt = y^2 - t
ode = Eq(Derivative(y, t), y^2 - t)

# Find the analytical solution
solution = dsolve(ode)

# Display the solution
println(solution)

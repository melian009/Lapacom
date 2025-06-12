using Pkg
using Symbolics
using LinearAlgebra
using DifferentialEquations
using DiffEqParamEstim
using Optim
using Plots


#=Test 1

function f(x)
    F = zero.(x)
    F[1] = x[1]^2 + x[3]
    F[2] = x[1] + x[2]
    F[3] = x[2]^2 + x[3]^2

    return F
end

x0 = [1,2,3]
J0 = jacobian(f, x0)

#Output
#3×3 Matrix{Int64}:
# 2  0  1
# 1  1  0
# 0  4  6

#Test 2
Symbolics.jacobian([x + x*y, x^2 + y],[x, y])
=#

using Symbolics

@variables N S X H r R K d gamma Smax

# Definición de las ecuaciones diferenciales
dNdt = X * r * N * R * (K - N/K) - d * N - H * N
dSdt = gamma * (S - (S^2) / (Smax - Smax * (1 - X) * H))

# Encontrar los puntos de equilibrio
equilibrium_points = solve([dNdt, dSdt], [N, S])

# Calcular la matriz Jacobiana
jac = jacobian([dNdt, dSdt], [N, S])

# Calcular el determinante de la matriz Jacobiana
det_jac = det(jac)

# Despejar la función H a partir del determinante de la matriz Jacobiana
H_eq = solve(det_jac)

# Imprimir los resultados
println("Puntos de Equilibrio:")
for eq_point in equilibrium_points
    println(eq_point)
end

println("Matriz Jacobiana:")
println(jac)

println("Determinante de la matriz Jacobiana:")
println(det_jac)

println("Función H en función del determinante de la matriz Jacobiana:")
println(H_eq)

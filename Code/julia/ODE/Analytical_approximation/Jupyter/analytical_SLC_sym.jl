#https://cooperrc.github.io/Julia-learning/day_07.html

using Plots, DifferentialEquations, LaTeXStrings, Symbolics
using ForwardDiff
import ForwardDiff.jacobian

@variables x_1, x_2

f = [x_1^2 + x_2, x_2^2 + x_1]

f_x= build_function(f, x_1, x_2, expression = Val{false})

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


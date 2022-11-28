#https://cooperrc.github.io/Julia-learning/day_07.html

using Plots, DifferentialEquations, LaTeXStrings, Symbolics
using ForwardDiff
import ForwardDiff.jacobian

@variables x_1, x_2

f = [x_1^2 + x_2, x_2^2 + x_1]

f_x= build_function(f, x_1, x_2, expression = Val{false})

out = f_x[1](1.0, 2.0)

Plots.theme(:cooper)
x1_array = range(0, 1, length = 100)
x2_array = range(0, 1, length = 100)
Farray1 = [f_x[1](x1i, x2i)[1] for x1i in x1_array for x2i in x2_array]
Farray2 = [f_x[1](x1i, x2i)[2] for x1i in x1_array for x2i in x2_array]

l = @layout [a; b]
p = contourf(x1_array, x2_array, Farray1, layout = l)
contourf!(p[2], x1_array, x2_array, Farray2)


#Evaluate the Jacobian
DJ = Symbolics.jacobian(f, [x_1, x_2])
DJ_xy = build_function(DJ, [x_1, x_2], expression = Val{false})
DJ_xy[1]([3, 2])

#Check out output
#2×2 Matrix{Int64}:
# 6  1
# 1  4

#Success!






julia> @variables y1 y2 r K d a E g
8-element Vector{Num}:
 y1
 y2
  r
  K
  d
  a
   E
  g

julia> J = Symbolics.jacobian([(r * y2 * ((K - y2) / K) - (d * y1) - (g * y1)), (g * y1) - (a * y2) - (E * y2)],[y1, y2])
2×2 Matrix{Num}:
 -d - g  (r*(K - y2) - r*y2) / K
      g                        -E - a

julia> det(J)
(-E - a)*(-d - g) + (-g*(r*(K - y2) - r*y2)) / K

julia> A=det(J) 

julia> M = Symbolics.simplify(A)
(E*K*d1 + K*d1*d2 + E*K*g + K*d2*g + 2g*r*x2 - K*g*r) / K


# https://www.juliapackages.com/p/symbolics
# Despejar x2 del determinante, obtener la x1 de las escuaciones y hacer las simulaciones para todo el rango de E



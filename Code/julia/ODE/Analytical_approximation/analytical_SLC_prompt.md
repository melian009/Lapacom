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

julia> 

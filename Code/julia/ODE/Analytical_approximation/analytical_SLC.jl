using Pkg
Pkg.activate(".")
using ForwardDiff
using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using GlobalSensitivity
using Statistics
using DiffEqParamEstim
using Optim
using Plots
using Plots.PlotMeasures
using Symbolics
import ForwardDiff.jacobian
using NonlinearSolve


@variables Na Sa r X R d K H g Smax
J_SLC = Symbolics.jacobian([(X * r * R * Na)*((K - Na)/K) - (d * Na) - H*(1 - X)*Na, (g * Sa) * (1 - (Sa/(Smax * (1 - H * (1-X)))))], [Na, Sa])
#OUTPUT
# [Na, Sa])
#2×2 Matrix{Num}:
# (R*X*r*(K - Na) - Na*R*X*r) / K - d - H*(1 - X)                                                                          0
#                                               0  g*(1 + (-Sa) / (Smax*(1 - H*(1 - X)))) + (-Sa*g) / (Smax*(1 - H*(1 - X)))

@variables Na Nb Sa Sb X Ra Rb ra rb da db K H ga gb Smaxa Smaxb c
J_SLC = Symbolics.jacobian([(ra * Ra * Na)*((K - Na)/K) - (da * Na) - H*(1 - X)*Na - c*Na*Nb, (rb * Rb * Nb)*((K - Nb)/K) - (db * Nb) - H*(1 - X)*Nb - c*Na*Nb,  (ga * Sa) * (1 - (Sa/(Smaxa - Smaxa * H * (1 - X)))), (gb* Sb) * (1 - (Sb/(Smaxb - Smaxb * H * (1 - X))))], [Na, Nb, Sa, Sb])

@variables Na Nb Sa Sb X R r d K H g Smax c
J_SLC = Symbolics.jacobian([(r * R * Na)*((K - Na)/K) - (d * Na) - H*(1 - X)*Na - c*Na*Nb, (r * R * Nb)*((K - Nb)/K) - (d * Nb) - H*(1 - X)*Nb - c*Na*Nb,  (g * Sa) * (1 - (Sa/(Smax - Smax * H * (1 - X)))), (g* Sb) * (1 - (Sb/(Smax - Smax * H * (1 - X))))], [Na, Nb, Sa, Sb])

det_jac = det(J_SLC)

X_SLC = Symbolics.jacobian([-d + ((K - Na)*R*r - Na*R*r) / K - H*(1 - X) - Nb*c,-Na*c,-Nb*c,-d + ((K - Nb)*R*r - Nb*R*r) / K - H*(1 - X) - Na*c], [Na, Nb])

#OUTPUT
#-da + ((K - Na)*Ra*ra - Na*Ra*ra) / K - H*(1 - X) - Nb*c, -Na*c, 0, 0
#-Nb*c    ,-db + ((K - Nb)*Rb*rb - Nb*Rb*rb) / K - H*(1 - X) - Na*c,  0, 0
#0,0,(-Sa*g) / (Smaxa - H*Smaxa*(1 - X)) + g*(1 + (-Sa) / (Smaxa - H*Smaxa*(1 - X))),0
#0,0,0,(-Sb*gb) / (Smaxb - H*Smaxb*(1 - X)) + gb*(1 + (-Sb) / (Smaxb - H*Smaxb*(1 - X)))

#Overleaf
#\begin{eqnarray}\label{eq:1}
#\frac{dN_{a_i}}{dt} &=&  r_{i} R_{i} N_{a_{i}} (K-N_{a_i}/K) - d_{i} N_{a_{i}} - \mathcal{H}_{i} (1 - \mathcal{X})N_{a_{i}} - c_{ij} N_{a_i} N_{a_j} \nonumber 
#\\
#\frac{dS_{a_{i}}}{dt} &=& \gamma_{i} S_{a_i} \Big(1-\frac{S_{a_i}}{S_{max_i}-S_{max_i} \mathcal{H}_{i} (1-\mathcal{X})} \Big),
#\end{eqnarray}


det_jac = det(J_SLC)

#OUTPUT
#(g*(1 + (-Sa) / (Smax*(1 - H*(1 - X)))) + (-Sa*g) / (Smax*(1 - H*(1 - X))))*((R*X*r*(K - Na) - Na*R*X*r) / K - d - H*(1 - X))

#OUTPUT
(-da + ((K - Na)*Ra*ra - Na*Ra*ra) / K - H*(1 - X) - Nb*c)*(-db + ((K - Nb)*Rb*rb - Nb*Rb*rb) / K - H*(1 - X) - Na*c)*((-Sa*g) / (Smaxa - H*Smaxa*(1 - X)) + g*(1 + (-Sa) / (Smaxa - H*Smaxa*(1 - X))))*((-Sb*gb) / (Smaxb - H*Smaxb*(1 - X)) + gb*(1 + (-Sb) / (Smaxb - H*Smaxb*(1 - X)))) - Na*Nb*(c^2)*((-Sa*g) / (Smaxa - H*Smaxa*(1 - X)) + g*(1 + (-Sa) / (Smaxa - H*Smaxa*(1 - X))))*((-Sb*gb) / (Smaxb - H*Smaxb*(1 - X)) + gb*(1 + (-Sb) / (Smaxb - H*Smaxb*(1 - X))))

   
M = Symbolics.simplify(det_jac)

#OUTPUT
#((Smax*g + H*Smax*X*g - 2.0Sa*g - H*Smax*g)*(H*K*X + K*R*X*r - H*K - K*d - 2Na*R*X*r)) / (K*Smax*(1 + H*X - H))

Symbolics.solve_for([(X * r * R * Na)*((K - Na)/K) - d * Na - H*(1 - X)*Na, g * Sa * (1 - (Sa/(Smax * (1 - H * (1-X)))))], [Na,Sa]; simplify=true)

#OUTPUT

#EXAMPLE
#https://github.com/JuliaSymbolics/Symbolics.jl/issues/381
#@variables r y x₂ h₁ g₁ x₁ g₂ h₂ g₄ g₃ f;
#Symbolics.solve_for([(r - y - x₂ * h₁) * g₁ ~ x₁,(x₁ - h₂ * y) * g₂ ~ x₂,x₂ * g₃ + x₁ * g₄ + f ~ y], [x₁,x₂,y]; simplify=false)


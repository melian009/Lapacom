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

#Complex model
@variables Na Nb Sa Sb X Ra Rb ra rb da db K H ga gb Smaxa Smaxb c
J_SLC = Symbolics.jacobian([(ra * Ra * Na)*((1 - Na)/K) - (da * Na) - H*(1 - X)*Na - c*Na*Nb, (rb * Rb * Nb)*((1 - Nb)/K) - (db * Nb) - H*(1 - X)*Nb - c*Na*Nb,  (ga * Sa) * (1 - (Sa/(Smaxa - Smaxa * H * (1 - X)))), (gb* Sb) * (1 - (Sb/(Smaxb - Smaxb * H * (1 - X))))], [Na, Nb, Sa, Sb])

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

M = Symbolics.simplify(det_jac)
#too complexxxxxxxxx

#Simplified model
@variables Na Nb Sa Sb X R r d K H g Smax c
X_SLC = Symbolics.jacobian([(r * R * Na)*((1 - Na)/K) - (d * Na) - H*(1 - X)*Na - c*Na*Nb, (r * R * Nb)*((1 - Nb)/K) - (d * Nb) - H*(1 - X)*Nb - c*Na*Nb,  (g * Sa) * (1 - (Sa/(Smax - Smax * H * (1 - X)))), (g* Sb) * (1 - (Sb/(Smax - Smax * H * (1 - X))))], [Na, Nb, Sa, Sb])

#OUTPUT
#-d + ((K - Na)*R*r - Na*R*r) / K - H*(1 - X) – Nb*c	-Na*c	0	0
#-Nb*c	-d + ((K - Nb)*R*r - Nb*R*r) / K - H*(1 - X) - Na*c	0	0
#0	0	(-Sa*g) / (Smax - H*Smax*(1 - X)) + g*(1 + (-Sa) / (Smax - H*Smax*(1 - X)))	0
#0	0	0	(-Sb*g) / (Smax - H*Smax*(1 - X)) + g*(1 + (-Sb) / (Smax - H*Smax*(1 - X)))


#Block matrix = Na Nb vs. Sa Sb
X_SLC = Symbolics.jacobian([(r * R * Na)*((K - Na)/K) - (d * Na) - H*(1 - X)*Na - c*Na*Nb, (r * R * Nb)*((K - Nb)/K) - (d * Nb) - H*(1 - X)*Nb - c*Na*Nb], [Na, Nb])

#-d + ((K - Na)*R*r - Na*R*r) / K - H*(1 - X) - Nb*c,-Na*c
#-Nb*c,-d + ((K - Nb)*R*r - Nb*R*r) / K - H*(1 - X) - Na*c

det_jac = det(X_SLC)

#OUTPUT
#(-d + ((K - Nb)*R*r - Nb*R*r) / K - H*(1 - X) - Na*c)*(-d + ((K - Na)*R*r - Na*R*r) / K - H*(1 - X) - Nb*c) - Na*Nb*(c^2)


   
M = Symbolics.simplify(det_jac)
#((H^2)*(K^2) + 2H*(K^2)*d + (K^2)*(d^2) - 2(H^2)*(K^2)*X + H*(K^2)*Na*c + H*(K^2)*Nb*c - 2H*(K^2)*R*r - 2H*(K^2)*X*d + 2H*K*Na*R*r + 2H*K*Nb*R*r + (K^2)*Na*c*d + (K^2)*Nb*c*d - 2(K^2)*R*d*r + 2K*Na*R*d*r + 2K*Nb*R*d*r + (H^2)*(K^2)*(X^2) - H*(K^2)*Na*X*c - H*(K^2)*Nb*X*c + 2H*(K^2)*R*X*r - 2H*K*Na*R*X*r - 2H*K*Nb*R*X*r - (K^2)*Na*R*c*r - (K^2)*Nb*R*c*r + (K^2)*(R^2)*(r^2) + 2K*(Na^2)*R*c*r - 2K*Na*(R^2)*(r^2) + 2K*(Nb^2)*R*c*r - 2K*Nb*(R^2)*(r^2) + 4Na*Nb*(R^2)*(r^2)) / (K^2)

Symbolics.solve_for([((H^2)*(K^2) + 2H*(K^2)*d + (K^2)*(d^2) - 2(H^2)*(K^2)*X + H*(K^2)*Na*c + H*(K^2)*Nb*c - 2H*(K^2)*R*r - 2H*(K^2)*X*d + 2H*K*Na*R*r + 2H*K*Nb*R*r + (K^2)*Na*c*d + (K^2)*Nb*c*d - 2(K^2)*R*d*r + 2K*Na*R*d*r + 2K*Nb*R*d*r + (H^2)*(K^2)*(X^2) - H*(K^2)*Na*X*c - H*(K^2)*Nb*X*c + 2H*(K^2)*R*X*r - 2H*K*Na*R*X*r - 2H*K*Nb*R*X*r - (K^2)*Na*R*c*r - (K^2)*Nb*R*c*r + (K^2)*(R^2)*(r^2) + 2K*(Na^2)*R*c*r - 2K*Na*(R^2)*(r^2) + 2K*(Nb^2)*R*c*r - 2K*Nb*(R^2)*(r^2) + 4Na*Nb*(R^2)*(r^2)) / (K^2)
], [Na,Nb]; simplify=true)

#EXAMPLE
#https://github.com/JuliaSymbolics/Symbolics.jl/issues/381
#@variables r y x₂ h₁ g₁ x₁ g₂ h₂ g₄ g₃ f;
#Symbolics.solve_for([(r - y - x₂ * h₁) * g₁ ~ x₁,(x₁ - h₂ * y) * g₂ ~ x₂,x₂ * g₃ + x₁ * g₄ + f ~ y], [x₁,x₂,y]; simplify=false)


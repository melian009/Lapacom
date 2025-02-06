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

#EQ 5 ms 
#\begin{equation}\label{eq:1}
#\frac{dN_i}{dt} =  r_i R_i N_i (1-N_{A_i}/K) - \mu_i N_i - \mathcal{H}_i (\mathcal{X})N_i - N_i \sum_{j=i:2} c_{ij} N_j~,     
#\end{equation}

#Session
#https://chatgpt.com/share/67a4af63-be44-800a-a843-b4b6746e99a6
#https://discourse.julialang.org/t/modeling-and-solving-systems-of-nonlinear-dynamic-equations/108556/12

#-----------------------------------------------------------
# for Hi(X) == 0, reproduction is not occurring so species can be exploited
#The reproductive cycle is introduced as a binary factor, $\mathcal{X}$, which equals 1 when species reproduces and exploitation is forbidden and 0 when the species are exploited.

@variables N1 N2 r1 R1 r2 R2 K mu1 mu2 H1 H2 c11 c12 c21 c22

# Define the system of ODEs
f1 = (r1 * R1 * N1) * ((1 - N1) / K) - (mu1 * N1) - (H1 * N1) - (N1 * N2 * c11) - (N1 * N2 * c12) - (N1 * N2 * c21) - (N1 * N2 * c22)
f2 = (r2 * R2 * N2) * ((1 - N2) / K) - (mu2 * N2) - (H2 * N2) - (N1 * N2 * c11) - (N1 * N2 * c12) - (N1 * N2 * c21) - (N1 * N2 * c22)

# Compute the Jacobian matrix
J_SLC = Symbolics.jacobian([f1, f2], [N1, N2])

#OUTPUT 
#julia> J_SLC = Symbolics.jacobian([f1, f2], [N1, N2])
#2×2 Matrix{Num}:
# -H1 - mu1 + (-N1*R1*r1 + (1 - N1)*R1*r1) / K - N2*c11 - N2*c12 - N2*c21 - N2*c22                                                -N1*c11 - N1*c12 - N1*c21 - N1*c22
#                                               -N2*c11 - N2*c12 - N2*c21 - N2*c22  -H2 - mu2 + (-N2*R2*r2 + (1 - N2)*R2*r2) / K - N1*c11 - N1*c12 - N1*c21 - N1*c22

# Display the result
display(J_SLC)

#julia> det_jac = det(J_SLC) 
#-(-N1*c11 - N1*c12 - N1*c21 - N1*c22)*(-N2*c11 - N2*c12 - N2*c21 - N2*c22) + (-H2 - mu2 + (-N2*R2*r2 + (1 - N2)*R2*r2) / K - N1*c11 - N1*c12 - N1*c21 - N1*c22)*(-H1 - mu1 + (-N1*R1*r1 + #(1 - N1)*R1*r1) / K - N2*c11 - N2*c12 - N2*c21 - N2*c22)

#To find the steady states of the system, we need to solve for N1​ and N2​ where both f1 and f2​ are equal to zero:
#f1(N1,N2)=0
#f2​(N1​,N2​)=0


#-------------------------------------------------------

#---------------------------------------------------------------
# for Hi(X) == 1, reproduction is occurring so species can not be exploited
#The reproductive cycle is introduced as a binary factor, $\mathcal{X}$, which equals 1 when species reproduces and exploitation is forbidden and 0 when the species are exploited. 
@variables N1 N2 r1 R1 r2 R2 K mu1 mu2 c11 c12 c21 c22

# Define the system of ODEs
f1 = (r1 * R1 * N1) * ((1 - N1) / K) - (mu1 * N1) - (N1 * N2 * c11) - (N1 * N2 * c12) - (N1 * N2 * c21) - (N1 * N2 * c22)
f2 = (r2 * R2 * N2) * ((1 - N2) / K) - (mu2 * N2) - (N1 * N2 * c11) - (N1 * N2 * c12) - (N1 * N2 * c21) - (N1 * N2 * c22)

# Compute the Jacobian matrix
J_SLC = Symbolics.jacobian([f1, f2], [N1, N2])

# Display the result
display(J_SLC)

#OUTPUT
#julia> J_SLC = Symbolics.jacobian([f1, f2], [N1, N2])
#2×2 Matrix{Num}:
# -mu1 + (-N1*R1*r1 + (1 - N1)*R1*r1) / K - N2*c11 - N2*c12 - N2*c21 - N2*c22                                           -N1*c11 - N1*c12 - N1*c21 - N1*c22
#                                          -N2*c11 - N2*c12 - N2*c21 - N2*c22  -mu2 + (-N2*R2*r2 + (1 - N2)*R2*r2) / K - N1*c11 - N1*c12 - N1*c21 - N1*c22

# det_jac = det(J_SLC) 
#julia> det_jac = det(J_SLC) 
#-(-N1*c11 - N1*c12 - N1*c21 - N1*c22)*(-N2*c11 - N2*c12 - N2*c21 - N2*c22) + (-mu2 + (-N2*R2*r2 + (1 - N2)*R2*r2) / K - N1*c11 - N1*c12 - N1*c21 - N1*c22)*(-mu1 + (-N1*R1*r1 + (1 - N1)*R1*r1) / K - N2*c11 - N2*c12 - N2*c21 - N2*c22)

#To find the steady states of the system, we need to solve for N1​ and N2​ where both f1 and f2​ are equal to zero:
#f1(N1,N2)=0
#f2​(N1​,N2​)=0

#--------------------------------------------------------------

#TODO
M = Symbolics.simplify(det_jac)
Symbolics.solve_for([det_jac], [N1,N2]; simplify=true)
symbolic_linear_solve([eq],[a])



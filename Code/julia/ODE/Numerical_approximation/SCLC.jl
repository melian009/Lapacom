#Packages
using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using GlobalSensitivity
using CairoMakie
using Statistics
using DataFrames
using CSV
using DiffEqParamEstim
using Optim
using Symbolics
using Plots
import ForwardDiff.jacobian

#=
Formulation of the simple life cicle for one site:

 Population dynamic:

    Simple Life cicle with only one class of population (only one class): Na + Sa

    dNa/dt = X * r * Ne * R * (K - Na/K) - da * Na - (1 - X) * H * Na - da * Na
    +
    dSa/dt = gamma * Sa * (1 - (Sa / (Smax - Smax * H * X)))

    Simple Life Cicle formulation that presents two clases of adult abundances:
    inmature and mature adults: Na, Nm + Sa

    dNa/dt = X * r * Ne * R * (K - Na/K) - gAM * Na - da * Na
    dNm/dt = gAM * Na * (K - (Na + Nm)/K) - (1 - X) * H * Nm - da * Nm
    +
    dSa/dt = gamma * Sa * (1 - (Sa / (Smax - Smax * H * X)))

Formulation for complex life cicles for one site: different numbers of clases.

  Population dinamic:

    Complex Life Cicle with 2 clases of population: Ne, Na + Sa
    
    dNe/dt = X * r * Na * R - gEA * Ne - de * Ne
    dNa/dt = gEA * Ne * (K - Na/K) - (1 - X) * H * Na - da * Na
    +
    dSa/dt = gamma * Sa * (1 - (Sa / (Smax - Smax * H * X)))

    Complex Life Cicle with 3 clases of population: Ne, Na, Nm + Sa

    dNe/dt = X * r * Na * R - gEA * Ne - de * Ne
    dNa/dt = gEA * Ne * (K - Na/K) - gAM * Na - da * Na
    dNm/dt = gAM * Na * (K - Nm/K) - (1 - X) * H * Nm - da * Nm
    +
    dSa/dt = gamma * Sa * (1 - (Sa / (Smax - Smax * H * X)))

    Complex Life Cicle with 5 clases of population: Ne, Nt, Nv, Nj, Na + Sa

    dNe/dt = X * r * Na * R - gEA * Ne - de * Ne
    dNt/dt = gET * Ne - gTV * Nt - dt * dNt
    dNv/dt = gTV * Nt * (K - Nv/K) - gVJ * Nv - dv * Nv
    dNj/dt = gVJ * Nv * (K - Nj/K) - gJA * Nj - dj * Nj
    dNa/dt = gJA * Nj * (K - Na/K) - (1 - X) * H * Na - da * Na
    +
    dSa/dt = gamma * Sa * (1 - ( Sa / (Smax - Smax * H * X)))

    Complex Life Cicle with 6 clases of population: Ne, Nt, Nv, Nj, Na, Nm + Sa

    dNe/dt = X * r * Na * R - gEA * Ne - de * Ne
    dNt/dt = gET * Ne - gTV * Nt - dt * dNt
    dNv/dt = gTV * Nt * (K - Nv/K) - gVJ * Nv - dv * Nv
    dNj/dt = gVJ * Nv * (K - Nj/K) - gJA * Nj - dj * Nj
    dNa/dt = gJA * Nj * (K - Na/K) - gAM * Na - da * Na
    dNm/dt = gAM * Na * (K - (Na + Nm)/K) - (1 - X) * H * Na - da * Na
    +
    dSa/dt = gamma * Sa * (1 - ( Sa / (Smax - Smax * H * X)))

 where i = nº of species == ["Patella ordinaria", "Patella aspera"]

Parameters and variables:
 - Ne = eggs abundance
 - Nt = trocophore abuncance
 - Nv = veliger abuncance
 - Nj = juvenile abuncance
 - Na = adults abuncance (when Nm is in the equations, Na means non matured adults)
 - Nm = matured adult abundaces
 - Sa = adults size (Average sizes Before MPS+FULL = [46.44,44.45])
 - r = population growth rate [9.17,5.03]
 - R = reproductive capacity
 - K = carrying capacity (k =1e^4)
 - X = Reproductive period [1,0] 
 - (1-X) = Exploitation periosd [1,0]
 - H = Exploitation rate (H = [0.639,0.57])
 - gEA = instant conversion rate of life stages (EA = Eggs to Adults) (gEA = 0.006)
 - de = natural mortality rate or death rate for eggs (de = [de_po,de_pa]) # Note: this value needs to be defined.
 - da = natural mortality rate or death rate for adults (da = [0.55,0.59]) # Note: empirical estimated values
 - Smax = maximum adult size estimated (56.0mm) # Note: Empirical value from the sample analized
 - gamma = adult growth rate (gamma=[0.32,0.36] year^{-1})

Empirical estimated mortalities (Z,d,F) and exploitation rates (E):

 {Patella ordinaria} (Henriques, 2012)

 - Natural Mortality rate (d) was 0.55 year^{-1}
 - Fishing mortality rate (F) was 1.24 year^{−1} 
 - Total mortality rate (Z=d+F) was 1.79 year^{−1}; 
 - Exploitation rate (E=F/Z) was 0.693. 

 {Patella aspera} (Sousa, 2017)

 - Natural mortality (d) was 0.59 year^{-1};
 - Fishing mortality (F) was 0.79 year^{-1};
 - Total mortality (Z=d+F) was 1.38 year^{-1};
 - Exploitation rate (E=F/Z) was 0.57. 

 For numerical simulation use these rates (d and E).

Average sizes before and after marine protected area implementations

 Before (1996-2006): FULL ACCESS.
 Patella apera = 43.53mm
 Patella ordinaria = 46.26mm

 After  (2007-2017): MPA+FULL
 Patella aspera = 44.45mm
 Patella ordinaria = 46.44mm

 Only MPS (2007-2017)
 Patella aspera = 50.61mm
 Patella ordinaria = 49.25mm

 Only Full Access (2007-2017) 

 Patella aspera = 43.41mm
 Patella ordinaria = 45.72mm
=#

# Full access scenario

function SLC!(du, u, p, t)
   Na, Sa = u
   r, R, K, H, X, da, Smax = p
  du[1] = dNa = X(t) * r * Na * R * (K - Na/K) - (1 - X(t)) * H(i) * Na - (da[i] * Na) 
  du[2] = dSa = gamma * Sa * (1 - Sa / (Smax - (1 - H[i])))
end

function aSLC!(du, u, p, t)
  Na,Nm, Sa = u
  r, R, K, H, X, da, Smax = p
 du[1] = dNa = X(t) * r * Na * R * (K - Na/K) - (da * Na) 
 du[2] = dNm = gAM * Na * (K - (Na + Nm)/K) - (1 - X(t)) * H * Nm - (da * Nm) 
 du[3] = dSa = gamma * Sa * (1 - Sa / (Smax - (1 - H)))
end

function SCLC!(du, u, p, t)
  Ne, Na, Sa = u
   r, R, K, H, X, gEA, de, da, Smax, gamma = p
  du[1] = dNe = X(t) * r * Na * R - (de * Ne) - (gEA * Ne)
  du[2] = dNa = gEA * Ne * (K - Na/K) - (da * Na) - ((1 - X(t))* H * Na)
  du[3] = dSa = gamma * Sa * (1 - Sa / (Smax - (1 * (1 - X(t)) * H)))
end

function aSCLC!(du, u, p, t)
  Ne, Na, Sa = u
   r, R, K, H, X, gEA,gAM, de, da, Smax, gamma = p
  du[1] = dNe = X(t) * r * Na * R - (de * Ne) - (gEA * Ne)
  du[2] = dNa = gEA * Ne * (K - Ne/K) - gAM * Na - da * Na 
  du[3] = dNm = gAM * Na * (K - (Na + Nm)/K) - (1 - X(t)) * H * Nm - (da * Nm) 
  du[4] = dSa = gamma * Sa * (1 - Sa / (Smax - (1 * (1 - X(t)) * H)))
end

function CLC!(du, u, p, t)
  Ne, Nt, Nv, Nj, Na, Sa = u
   r,R,K, H, X, gET, gTV, gVJ, gJA, de, da, Smax, gamma = p
  du[1] = dNe = X(t) * r * Na * R- de * Ne - gEA * Ne
  du[2] = dNt = gET * Ne - gTV*Nt - dt * Nt
  du[3] = dNv = gTV * Nt * (K - Nt/K) - gVJ * Nv - dv * Nv
  du[4] = dNj = gVJ * Nv * (K - Nj/K) - gJA * Nj - da * Na
  du[5] = dNa = gJA * Nj * (K - Na / K) - da * Na - (1 - X(t)) * H * Na
  du[6] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - X(t))))
end

function aCLC!(du, u, p, t)
  Ne, Nt, Nv, Nj, Na, Sa = u
   r, R, K, H, X, gET, gTV, gVJ, gJA, de, da, Smax, gamma = p
  du[1] = dNe = X(t) * r * Na * R(Sa,Smax) - de * Ne - gEA * Ne
  du[2] = dNt = gET * Ne - gTV*Nt - dt * Nt
  du[3] = dNv = gTV * Nt * (K - Nt/K) - gVJ * Nv - dv * Nv
  du[4] = dNj = gVJ * Nv * (K - Nj/K) - gJA * Nj - da * Na
  du[5] = dNa = gJA * Nv * (K - Na/K) - gAM * Na - da * Na 
  du[6] = dNm = gAM * Na * (K - (Na + Nm)/K) - (1 - X(t)) * H * Nm - (da * Nm) 
  du[7] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - X(t))))
end

function R(Sa,Smax)
  Smat = 1.34 * Sa - 28.06
  R = min(max(0.5(1 + (Sa - Smaturity)/(Smax - Smaturity)), 0), 1)
end
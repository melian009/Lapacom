
#Packages

using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics


# Simple life cycle on one site.


#=Check analytical_SLC_without_S.jl 

#Exp_lim = 1                 # Exploitation max limit 
#m=0.05                      # Interval of exploitation values 
#for H = 0:m:Exp_lim         #exploitation gradient
#end
function MSLC!(du, u, t, p)
    Na, Sa,= u
    X, re, K, H, R, d, c12, xj, Smax, g = p
    
  
    du[1] = dNa = (X(t) * re * R * Na)*((K - Na)/K) - (d * Na) - H*(1 - X(t)) * Na - c12 * Na * x2
    du[2] = dSa = (g * Sa) * (1 - (Sa/(Smax * (1 - H * (1 - X(t))))))
end


function X_(t)
  phi(t) = 2*pi*(t-0)/0.42
  X(t) = 0.5*(1*-tanh(0.5-sen(phi(t))))
  return X(t)
end



u0 = [1000.0, 25.0]  # Initial conditions
tspan = (1.0, 365*2)  # Time span for the simulation (from t=0 to t=1000)
p_1=[X_,0.32, 640000, 0.639, 1, 0.55, 0.5, 10000, 56, 0.6]

prob = ODEProblem(MSLC!, u0, tspan, p_1)
sol = solve(prob, Tsit5()) # "Error: BoundsError: attempt to access Float64 at index [2]"

plot(sol, xlabel="Time", ylabel="State Variables", label=["x1" "x2" "y1" "y2"])


u0 = [2500.0,25.0]  # Initial conditions for [x, y]
tspan = (0.0, 1000)  # Time span for the simulation (from t=0 to t=1000)
period = zeros(Float64,size(1:365*2))
# X, R, r, K, d, H, c12, Smax
#Plot de t_ vs X_ para ver los periodos de reproducción intermitentes en el tiempo.


function X_(t,t_0,k)
  phi(t) = 2*pi*(t-t_0)/(0.42*365)
  X(t) = 0.5(1-sin(k-phi(t)))
  return X(t)
end

function X(t)
  if (t % 365) / 365 >= 0.42
    return 1.0 # Reproductive Cycle
  else
    return 0.0 # Exploitation Cycle
  end
end

period = zeros(Float64,size(1:365.14*5))
periodX = zeros(Float64,size(1:365.14*5))
c=1
for t_ in 1:(365.14*5)
period[c] = X_(t_,365*0.42,0.1)
periodX[c] = X(t_)
c=c+1
end
period 



plot(period, label="Reproductive cycle by gradient")
plot!(periodX, label="Reproductive cicle by factor")
xlabel!("Time (days)")
ylabel!("Reproductive Cycle")

savefig!("X_cycle.png", dpi = 300)
=#






# Non-trivial solution for Na and Sa on a Single Site.
#Parameters
avg_oocytes = mean([92098, 804183]) # This is the actual mean. 
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r_= reggs
K_ = 64000    # Carrying capacity
d1_ = 0.590
H1_ = 0.639
#c12_ = 0.5 #competition term species 2 on 1
Smax_ = 56.0             # Maximum size for adults
Naj = 2500
avg_size = 33.4




#Vectors for ploting and simulations
n=2    #Number of years in the simulation
t_span = length(zeros(Float64,size(1:365.14*n)))
h_span = length(zeros(Float64, size(0:0.1:1)))
span = ones(Float64,size(1:365.14*n))
Kspan = ones(Float64,size(1:365.14*n))*K_ 
Sm_span = ones(Float64,size(1:365.14*n))*Smax_/2   # Linea de Smax/2 

H_r = range(0, 1, length=h_span)
H_span = ones(Float64,h_span)
cij_span = ones(Float64,h_span)

zero_line = zeros(Float64,h_span)

for i in 1:length(H_span)
    H_span[i] = H_r[i]
    cij_span[i] = H_r[i]
end

H_span
cij_span

#Outputs of simulations
Nai_h = zeros(t_span,h_span)
Sai_h = zeros(t_span,h_span)
periodX = zeros(t_span)



for i in 1:length(cij_span)
  cij_ = cij_span[5]
  H_i = H_span[i]
  c=1 

  Nai1 = zeros(Float64,size(1:365.14*n))
  Sai1 = zeros(Float64,size(1:365.14*n))
  periodX = zeros(t_span)

  for t_ in 1:(365.14*n)
    t_0 = (365.14*0.42)
    k=0.1
      phi(t_) = 2*pi*(t_ - t_0)/(0.42 * 365)
      periodX[c] = tanh(1 - sin(k-phi(t_)))
      Smat = 1.34 * (avg_size) - 28.06
      R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax_ - Smat)), 0.0), 1.0)
      #Non trivial solution expresions
        Nai1[c] = -(H_i * (1 - periodX[c]) + cij_span[i] * Naj - periodX[c] * r_ * R_ + d1_)* K_ / (periodX[c] * r_ * R_ * 2) 
        Sai1[c] = 1/2*(Smax_ * (1 - H_i * (1 - periodX[c])))
    
    avg_size = mean(Sai1)  
    c=c+1
    end
    Nai_h[:,i] = Nai1
    Sai_h[:,i] = Sai1

end

#Sa vs t (by H)
plot(Nai_h, label=["H=0" "H=0.1" "H=0.2" "H=0.3" "H=0.4" "H=0.5" "H=0.6" "H=0.7" "H=0.8" "H=0.9" "H=1"], ylim=(0,K_/2*1.2), legend=:right)
plot!(Kspan/2,label=false, color=:red)
annotate!(120, K_/2*1.03, text("K/2", :red, :center, 7))
ylabel!("Adult abundances (nº individuals)")
xlabel!("Time (days)")
plot!(periodX*K_/2, label="X(t)",color=:blue, style = :dash)
plot!(span*K_/2, c=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, K_/2+1000, text("X(t)=1", :blue, :center, 8))
annotate!(50, -1000, text("X(t)=0", :blue, :center, 8))
title!("SLC Adult dinamics in differtent grades of H")


#Na vs t (by H)
plot(Sai_h, label=["H=0" "H=0.1" "H=0.2" "H=0.3" "H=0.4" "H=0.5" "H=0.6" "H=0.7" "H=0.8" "H=0.9" "H=1"], legend=:right)
plot!(Sm_span,color=:red,label = false)
annotate!(-55, Smax_/2, text("Smax/2", :red, :center, 8))
plot!(periodX*Smax_/2, label="X(t)",color=:blue, style = :dash)
plot!(span*Smax_/2, color=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, Smax_/2+1, text("X(t)=1", :blue, :center, 8))
annotate!(50, 1, text("X(t)=0", :blue, :center, 8))
xlabel!("Time (days)")
ylabel!("Adult Size (mm)")
title!("Adult size in different grades of H")



#Na and Sa vs H

Sa_H = ones(Float64,h_span)
Na_H = ones(Float64,h_span)

for j in 1:h_span
Sa_H[j] = minimum(Sai_h[:,j])
Na_H[j] = maximum(Nai_h[:,j])
end
Sa_H[1]
Smax_/2

#Sa vs H
plot(H_span,Sa_H, label=false)
xlabel!("Exploitation normalized rate (H)")
ylabel!("Asult Size (mm)")



#Na vs H
plot(H_span,Na_H,label=false)
plot!(H_span,zero_line,label=false)
xlabel!("Exploitation normalizedrate (H)")
ylabel!("Adult abundances (nº individuals)")






n=10    #Number of years in the simulation
t_span = length(zeros(Float64,size(1:365.14*n)))
h_span = length(zeros(Float64, size(0:0.1:1)))
span = ones(Float64,size(1:365.14*n))
Kspan = ones(Float64,size(1:365.14*n))*K_ 
Sm_span = ones(Float64,size(1:365.14*n))*Smax_/2   # Linea de Smax/2 


Nai_cij = zeros(t_span,h_span)
Sai_cij = zeros(t_span,h_span)
periodX = zeros(t_span)


H_r = range(0, 1, length=h_span)
H_span = ones(Float64,h_span)
cij_span = ones(Float64,h_span)

for i in 1:length(H_span)
  H_span[i] = H_r[i]
  cij_span[i] = H_r[i]
end




#cij gradient. 

for i in 1:length(H_span)
  cij_ = cij_span[i]
  H_i = H1_
  c=1 
  Nai1 = zeros(Float64,size(1:365.14*n))
  Sai1 = zeros(Float64,size(1:365.14*n))
  periodX = zeros(t_span)

  for t_ in 1:(365.14*n)
    t_0 = (365.14*0.42)
    k=0.1
      phi(t_) = 2*pi*(t_ - t_0)/(0.42 * 365)
      periodX[c] = tanh(1 - sin(k-phi(t_)))
      Smat = 1.34 * (avg_size) - 28.06
      R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax_ - Smat)), 0.0), 1.0)
      #Non trivial solution expresions
        Nai1[c] = -(H_i * (1 - periodX[c]) + cij_span[i] * Naj - periodX[c] * r_ * R_ + d1_)* K_ / (periodX[c] * r_ * R_ * 2) 
        Sai1[c] = 1/2*(Smax_ * (1 - H_i * (1 - periodX[c])))
    
    avg_size = mean(Sai1)  
    c=c+1
    end
    Nai_cij[:,i] = Nai1
    Sai_cij[:,i] = Sai1

end


#Na vs t (by cij)
plot(Nai_cij, label=["cij=0" "cij=0.1" "cij=0.2" "cij=0.3" "cij=0.4" "cij=0.5" "cij=0.6" "cij=0.7" "cij=0.8" "cij=0.9" "cij=1"], ylim=(0,K_/2*1.2), legend=:right)
plot!(Kspan/2,label=false, color=:red)
annotate!(120, K_/2*1.03, text("K/2", :red, :center, 7))
ylabel!("Abuncance population")
xlabel!("Time (days)")
title!("SLC Adult dinamics in differtent grades of cij (Naj=2500ind)")
plot!(periodX*K_/2, label="X(t)",color=:blue, style = :dash)
plot!(span*K_/2, c=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, K_/2+1000, text("X(t)=1", :blue, :center, 8))
annotate!(50, -1000, text("X(t)=0", :blue, :center, 8))



#Sa vs t (by cij)
plot(Sai_cij, label=["H=0" "H=0.1" "H=0.2" "H=0.3" "H=0.4" "H=0.5" "H=0.6" "H=0.7" "H=0.8" "H=0.9" "H=1"], legend=:right)
plot!(Sm_span,color=:red,label = false)
annotate!(100, Smax_/2+0.5, text("Smax/2", :red, :center, 8))
plot!(periodX*Smax_/2, label="X(t)",color=:blue, style = :dash)
plot!(span*(Smax_)/2 , color=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(120, 26, text("X(t)=1", :blue, :center, 8))
annotate!(50, 1, text("X(t)=0", :blue, :center, 8))
xlabel!("Time (days)")
ylabel!("Adult Size (mm)")
title!("SLC Adult size in differtent grades of cij")





Sa_c = ones(Float64,h_span)
Na_c = ones(Float64,h_span)

for j in 1:h_span
Sa_c[j] = maximum(Sai_cij[:,j])
Na_c[j] = maximum(Nai_cij[:,j])
end


#Sa vs cij
plot(H_span,Sa_c, label=false, ylim = (27.34,27.36)) 
xlabel!("Species interaction coeficent cij")
ylabel!("Asult Size (mm)")



#Na vs cij
plot(H_span,Na_c,label=false)
plot!(H_span,zero_line,label=false)
xlabel!("Species interaction coeficent cij")
ylabel!("Adult abundances (nº individuals)")


#Complex life cycle on a single site
#No da lo que debería, al principio tenía resultados coherentes, pero ahora me da error y no hace las simulaciones debidamente:
#Warning: dt(2.2737367544323206e-13) <= dtmin(2.2737367544323206e-13) at t=34.86795993721168, and step error estimate = 2.7127354535020856. Aborting. There is either an error in your model specification or the true solution is unstable.
#└ @ SciMLBase C:\Users\Usuario\.julia\packages\SciMLBase\McEqc\src\integrator_interface.jl:599
#retcode: DtLessThanMin

function CLC!(du, u, p, t)
  Ne, Nt, Nv, Nj, Na, Sa = u
  t_0, k, r, K, H , d, Smax, gamma = p
   
   
   
   #reproductive cycle
   phi(t) = 2*pi*(t - t_0)/(0.42 * 365)
   periodX(t) = tanh(1 - sin(k-phi(t)))

   # reproductive capacity
   avg_size = du[6]
   Smat = 1.34 * (avg_size) - 28.06
   R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax - Smat)), 0.0), 1.0)

  du[1] = dNe = periodX(t) * r[1] * Na * R_ - d[1] * Ne - r[2] * Ne
  du[2] = dNt = periodX(t)*r[2] * Ne - r[3] * Nt - d[2] * Nt
  du[3] = dNv = periodX(t)*r[3] * Nt * ((K - Nv)/ K) - r[4] * Nv - d[3] * Nv
  du[4] = dNj = periodX(t)*r[4] * Nv * ((K - Nj)/ K) - r[5] * Nj - d[4] * Nj
  du[5] = dNa = periodX(t)*r[5]* R_ * Nj * ((K - Na)/ K) - d[5] * Na - (1 - periodX(t)) * H * Na
  du[6] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - periodX(t))))

end

avg_oocytes = mean([92098, 804183]) # This is the actual mean. mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r_ = [reggs, 0.998611, 0.971057, 0.4820525, 0.00629]
# natural death rates per life stage.
d = [0.999/365.14, 0.585/365.14, 0.343/365.14, 0.201/365.14, 0.000322/365.14]  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = 0.00014898749263737575
t0_ = 365.14*0.42
k_= 0.1
cij_= 0.5
Naj = 2500

K_ = 64000.00    # Carrying capacity
H1_ = 0.1
Smax_ = 56.0             # Maximum size for adults

#         t_0, k,  r,  K,  H ,  d, Smax,  gamma, cij, xj
p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate] 



n=3   #Number of years in the simulation
tspan = (0,365.14*n)
U0_ = [10^7,0,0,0,10^4, 49.25]
prob_ = ODEProblem(CLC!, U0_, tspan, p_span)
solve_= solve(prob_, Tsit5())




t_l = length(zeros(Float64,size(1:length(solve_.u))))
n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))
vars = zeros(t_l,n_v)
time = zeros(t_l,1)
for j in 1:length(solve_.u)
  for i in 1:6
    vars[j,i] = solve_.u[j][i]
    time[j] = solve_.t[j]
   end 
end

Ne = vars[:,1]
Nt = vars[:,2]
Nv = vars[:,3]
Nj = vars[:,4]
Na = vars[:,5]
Sa = vars[:,3]


time2 = time .+ 0.7



plot(time, Ne, label="Ne",xlim = (800,1000))#, ylim=(6.5*10^4,2.5*10^7))
plot!(time .- 0.7, Nt, label="Nt",xlim = (800,1000))#, ylim=(6.5*10^4,2.5*10^7))
plot!(time .- (0.7+1.3), Nv, label="Nv",xlim = (800,1000), ylim=(0,6.5*10^4))
plot!(time .- (0.7+1.3 + 7), Nj, label="Nj",xlim = (800,1000))#, ylim=(0,6.5*10^4))
plot!(time .+ 230 , Na, label="Na",xlim = (800,1000))
plot!(Kspan)

plot(solve_, vars=6, label="Sa", color=:black)


t
phi(t) = 2*pi*(t - t_0)/(0.42 * 365)
periodX(t) = tanh(1 - sin(k-phi(t)))

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
    Na,
     Sa,= u
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
n=1.5    #Number of years in the simulation
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
      phi(t_) = 2*pi*(t_ - t_0)/(365.25)
      periodX[c] = 1/2*(1+tanh(2*sin(phi(t_)) - k))
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
HX_NAt = plot(Nai_h, label=["H=0" "H=0.1" "H=0.2" "H=0.3" "H=0.4" "H=0.5" "H=0.6" "H=0.7" "H=0.8" "H=0.9" "H=1"], ylim=(0,K_/2*1.2), legend =:outerright)
plot!(Kspan/2,label=false, color=:red)
annotate!(50, K_/2*1.03, text("K/2", :red, :center, 7))
ylabel!("Adult abundance (nº individuals)")
xlabel!("Time (days)")
plot!(periodX*(K_/2), label="X(t)",color=:blue, style = :dash)
plot!(span*(31299.796819466865), c=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, 30000, text("X(t)=1", :blue, :center, 8))
annotate!(50, 3000, text("X(t)=0", :blue, :center, 8))
title!("Na* in differtent grades of H")
#png("HX_NAt.png")

display(HX_NAt)

#Na vs t (by H)

HX_SAt = plot(Sai_h, label=["H=0" "H=0.1" "H=0.2" "H=0.3" "H=0.4" "H=0.5" "H=0.6" "H=0.7" "H=0.8" "H=0.9" "H=1"], legend =:outerright)
plot!(Sm_span,color=:red,label = false)
annotate!(50, Smax_/2+1, text("Smax/2", :red, :center, 7))
plot!(periodX*Smax_/2, label="X(t)",color=:blue, style = :dash)
plot!(span*(27.387286092247592), color=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(250, 26, text("X(t)=1", :blue, :center, 8))
annotate!(250, 1, text("X(t)=0", :blue, :center, 8))
xlabel!("Time (days)")
ylabel!("Adult Size (mm)")
title!("Sa* in different grades of H")
#png("HX_SAt.png")

display(HX_SAt)



#Na and Sa vs H

Sa_H = ones(Float64,h_span)
Na_H = ones(Float64,h_span)

for j in 1:h_span
Sa_H[j] = minimum(Sai_h[:,j])
Na_H[j] = maximum(Nai_h[:,j])
end
Sa_H
Na_H
Smax_/2

#Sa vs H
H_Sa_ = plot(H_span,Sa_H, label=false)
xlabel!("Exploitation rate (H)")
ylabel!("Adult Size (mm)")
#png("H_Sa_")

display(H_Sa_)


#Na vs H
H_Na_ = plot(H_span,Na_H,label=false, ylims=(0,32000))
annotate!(0.83, 700, text("|H = 0.77 ", :blue, :center, 8))
xlabel!("Exploitation rate (H)")
ylabel!("Adult abundance (nº individuals)")
#png("H_Na_")

display(H_Na_)



n=1.5    #Number of years in the simulation
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
N_H_c_span = zeros(11,11)
S_H_c_span = zeros(11,11)

for i in 1:length(H_span)
  H_span[i] = H_r[i]
  cij_span[i] = H_r[i]
end


for j in 1:length(H_span)
    H_j = H_span[j]
    cij_i = cij_span[5]
    c=1
    n=1.5 #Years
    Nai1 = zeros(Float64,size(1:365.14*n))
    Sai1 = zeros(Float64,size(1:365.14*n))
    periodX = zeros(t_span)
    avg_size = 33.4
    for t_ in 1:(365.25*n)
      t_0 = (365.25*0.42)
      k=0.42
        phi(t_) = 2*pi*(t_ - t_0)/(365.25)
        periodX[c] = 1/2*(1+tanh(2*sin(phi(t_)) - k))
        Smat = 1.34 * (avg_size) - 28.06
        R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax_ - Smat)), 0.0), 1.0)
      #Non trivial solution expresions
        Nai1[c] = -(H_j * (1 - periodX[c]) + cij_i * Naj - periodX[c] * r_ * R_ + d1_)* K_ / (periodX[c] * r_ * R_ * 2) 
        Sai1[c] = 1/2*(Smax_ * (1 - H_j * (1 - periodX[c])))
    
        avg_size = mean(Sai1)  
        c=c+1
    end
    Nai_cij[:,j] = Nai1
    Sai_cij[:,j] = Sai1
  end 
  Sa_c = ones(Float64,h_span)
  Na_c = ones(Float64,h_span)
  for m in 1:h_span
  Sa_c[m] = minimum(Sai_cij[:,m])
  Na_c[m] = maximum(Nai_cij[:,m])
end


for i in 1:length(cij_span)
  cij_i = cij_span[i]
  
  N_H_c_span[:,i] = Na_c
  S_H_c_span[:,i] = Sa_c

end

N_H_c_span

#Na vs t (by cij)
cijX_NAt = plot(Nai_cij, label=["cij=0" "cij=0.1" "cij=0.2" "cij=0.3" "cij=0.4" "cij=0.5" "cij=0.6" "cij=0.7" "cij=0.8" "cij=0.9" "cij=1"], ylim=(0,K_/2*1.2), legend=:outerright)
plot!(Kspan/2,label=false, color=:red)
annotate!(0, K_/2*1.03, text("K/2", :red, :center, 7))
ylabel!("Abuncance population (Nº of individuals)")
xlabel!("Time (days)")
title!("Na* in differtent grades of cij")
plot!(periodX*(K_/2), label="X(t)",color=:blue, style = :dash)
plot!(span*(31299.796819466865), c=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, 30000, text("X(t)=1", :blue, :center, 8))
annotate!(50, 3000, text("X(t)=0", :blue, :center, 8))
#png("cijX_NAt")


#Sa vs t (by cij)
cijX_SAt = plot(Sai_cij, label=["cij=0" "cij=0.1" "cij=0.2" "cij=0.3" "cij=0.4" "cij=0.5" "cij=0.6" "cij=0.7" "cij=0.8" "cij=0.9" "cij=1"], legend=:outerright)
plot!((Sm_span),color=:red,label = false)
annotate!(30, (Smax_)/2+0.5, text("Smax/2", :red, :center, 8))
plot!(periodX*(Smax_+1)/2, label="X(t)",color=:blue, style = :dash)
plot!(span*(Smax_)/2 , color=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(250, 26, text("X(t)=1", :blue, :center, 8))
annotate!(250, 1, text("X(t)=0", :blue, :center, 8))
xlabel!("Time (days)")
ylabel!("Adult Size (mm)")
title!("Sa* in differtent grades of cij")
#png("cijX_SAt")




Sa_c = ones(Float64,h_span)
Na_c = ones(Float64,h_span)

for j in 1:h_span
Sa_c[j] = minimum(Sai_cij[:,j])
Na_c[j] = maximum(Nai_cij[:,j])
end


#Sa vs cij
cij_Sa = plot(H_span,Sa_c, label=false, ylims=(10.3,10.4))
annotate!(-0.07, 10.3723388021333, text("10.372 -", :blue, :center, 8))
xlabel!("Spatial competition coeficent cij")
ylabel!("Adult Size (mm)")
#png("cij_Sa")


#Na vs cij
cij_Na = plot(H_span,Na_c,label=false, ylims= (0,32000))
annotate!(0.835, 700, text("|cij = 0.77 ", :blue, :center, 8))
xlabel!("Spatial competition coeficent cij")
ylabel!("Adult abundances (nº individuals)")
#png("cij_Na")

#Complex life cycle on a single site
#No da lo que debería, al principio tenía resultados coherentes, pero ahora me da error y no hace las simulaciones debidamente:
#Warning: dt(2.2737367544323206e-13) <= dtmin(2.2737367544323206e-13) at t=34.86795993721168, and step error estimate = 2.7127354535020856. Aborting. There is either an error in your model specification or the true solution is unstable.
#└ @ SciMLBase C:\Users\Usuario\.julia\packages\SciMLBase\McEqc\src\integrator_interface.jl:599
#retcode: DtLessThanMin

function CLC!(du, u, p, t)
  Ne, Nt, Nv, Nj, Na, Sa = u
  t_0, k, r, K, H , d, Smax, gamma, cij, Naj = p
   
   
   
   #reproductive cycle
   phi(t) = 2*pi*(t - t_0)/(365)
   periodX(t) = 1/2*(1+tanh(2*sin(phi(t)) - k))

   # reproductive capacity
   avg_size = du[6]
   Smat = 1.34 * (avg_size) - 28.06
   R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax - Smat)), 0.0), 1.0)

  du[1] = dNe = periodX(t) * r[1] * Na * R_ - d[1] * Ne - r[2] * Ne
  du[2] = dNt = periodX(t)*r[2] * Ne - r[3] * Nt - d[2] * Nt
  du[3] = dNv = periodX(t)*r[3] * Nt*((K - Nv)/ K) - r[4] * Nv - d[3] * Nv
  du[4] = dNj = periodX(t)*r[4] * Nv * ((K - Nj)/ K) - r[5] * Nj - d[4] * Nj
  du[5] = dNa = periodX(t)*r[5]* R_ * Nj * ((K - Na)/ K) - d[5] * Na - (1 - periodX(t)) * H * Na - cij * Na * Naj
  du[6] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - periodX(t))))

end

avg_oocytes = mean([92098, 804183]) # This is the actual mean. mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r_ = [reggs, 0.998611, 0.971057, 0.4820525, 0.00629]
# natural death rates per life stage.
d = [0.999/365.25, 0.585/365.25, 0.343/365.25, 0.201/365.25, 0.000322/365.25]  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = 0.00014898749263737575
t0_ = 365.14*0.42
k_ = 0.42
K_ = 64000.00    # Carrying capacity
H1_ = 0
Smax_ = 56.0             # Maximum size for adults
cij_ = 0.5
Naj = 500
#         t_0, k,  r,  K,  H ,  d, Smax,  gamma, cij, xj


p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate, cij_, Naj] 
n=10   #Number of years in the simulation
tspan = (0,365.14*n)
U0_ = [10^7,0,0,0,10^4, 49.25]
prob_ = ODEProblem(CLC!, U0_, tspan, p_span)
solve_= solve(prob_, Tsit5())




t_l = length(zeros(Float64,size(1:length(solve_.u))))
n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))
vars = zeros(t_l,n_v)
days = ones(t_l,1)

for j in 1:length(solve_.u)
  for i in 1:6
    vars[j,i] = solve_.u[j][i]
    days[j] = solve_.t[j]
   end 
end

Ne = vars[:,1]
Nt = vars[:,2]
Nv = vars[:,3]
Nj = vars[:,4]
Na = vars[:,5]
Sa = vars[:,6]

plot!(days, Na, label="cij=0.5", legend=:outerright)






CLC_NeNt_ = plot(days.-(0.7 + 1.3 + 7 + 230), Ne, label="Eggs", xlim=(0,365.25*1.5), legend=:outerright, color=:red)
plot!(days.-(1.3 + 7 + 230), Nt, label="Trochophore", color=:blue)
xlabel!("Time (days)")
ylabel!("Abundance (nº individuals)")
png("CLC_NeNt_")


CLC_NAt0 = plot(days, Na, label="Adults", color = :black)
plot!(days.-(7 + 230), Nv, label="Veliger", legend=:outerright, color=:green)
plot!(days.-(230), Nj, label="Juvenile",ylims=(0,65000),xlims=(0,365.25*1.5), legend=:outerright,color=:brown)
xlabel!("Time (days)")
ylabel!("Abundance (nº individuals)")
png("CLC_NAt0")

CLC_NAt1 = plot(time.- (0.7 + 1.3 + 7 + 230), Ne, label="Ne", xlim=(150, 550), legend=:outerright, color=:red)
plot!(time.- (1.3 + 7 + 230), Nt, label="Nt", color=:blue)
xlabel!("time (days)")
ylabel!("Abundance (nº individuals)")
png("CLC_NAt1")

CLC_NAt2 = plot(time.-(1.3 + 7 + 230), Nt, label="Nt", xlim=(150,550),legend=:outerright, color=:blue)
plot!(time.-(7 + 230), Nv, label="Nv", legend=:outerright, color=:green)
xlabel!("time (days)")
ylabel!("Abundance (nº individuals)")
png("CLC_NAt2")

CLC_NAt3 = plot(time.-(7 + 230), Nv, label="Nv", xlim=(150,550),legend=:outerright, color=:green)
plot!(time.-(230), Nj, label="Nj",ylims=(0,65000), legend=:outerright,color=:brown)
xlabel!("time (days)")
ylabel!("Abundance (nº individuals)")
png("CLC_NAt3")

CLC_NAt4 = plot(time.-(230), Nj, label="Nj",ylims=(0,65000), xlim=(150,550), legend=:outerright,color=:brown)
plot!(time  , Na, label="Na", color = :black)
xlabel!("time (days)")
ylabel!("Abundance (nº individuals)")
png("CLC_NAt4")




CLC_SAt =  plot(time, Sa, label="Sa", color = :blue, legend=:outerright)
xlabel!("time (days)")
ylabel!("Adult size (mm)")
png("CLC_SAt")









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
n=10    #Number of years in the simulation
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
Nai_h_c = zeros(h_span,h_span)
Sai_h_c = zeros(h_span,h_span)


for m in 1:length(cij_span)
  for i in 1:length(H_span)
    cij_ = cij_span[m]
    H_i = H_span[i]
    c=1 

    Nai1 = zeros(Float64,size(1:365.14*n))
    Sai1 = zeros(Float64,size(1:365.14*n))
    periodX = zeros(t_span)

    for t_ in 1:(365.14*n)
      t_0 = (365.14*0.42)
      k=0.1
      phi(t_) = 2*pi*(t_ - t_0)/(365.25)
      periodX[c] = 1/2*(1+tanh(2*sin(phi(t_)) - k))
      Smat = 1.34 * (avg_size) - 28.06
      R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax_ - Smat)), 0.0), 1.0)
      #Non trivial solution expresions
        Nai1[c] = -(H_i * (1 - periodX[c]) + cij_ * Naj - periodX[c] * r_ * R_ + d1_)* K_ / (periodX[c] * r_ * R_ * 2) 
        Sai1[c] = 1/2*(Smax_ * (1 - H_i * (1 - periodX[c])))
    
      avg_size = mean(Sai1)  
      c=c+1
    end
    Nai_h[:,i] = Nai1
    Sai_h[:,i] = Sai1

  end
  #Na and Sa vs H

  Sa_H = ones(Float64,h_span)
  Na_H = ones(Float64,h_span)

  for j in 1:h_span
  Sa_H[j] = minimum(Sai_h[:,j])
  Na_H[j] = maximum(Nai_h[:,j])
  end
 Sai_h_c[:,m] = Sa_H
 Nai_h_c[:,m] = Na_H
end
 
real(log(Nai_h_c))

plot(H_span,real(log(Nai_h_c)), label=["cij=0" "cij=0.1" "cij=0.2" "cij=0.3" "cij=0.4" "cij=0.5" "cij=0.6" "cij=0.7" "cij=0.8" "cij=0.9" "cij=1"], legend=:outerright)
xlabel!("Exploitation rate (H)")
ylabel!("Adult abundance (nº individuals)")
png("H_decay_for_cij_gradient_N_log")

span_11 = ones(Float64,size(0:0.1:1))
plot(H_span, Nai_h_c,
      label=["cij=0" "cij=0.1" "cij=0.2" "cij=0.3" "cij=0.4" "cij=0.5" "cij=0.6" "cij=0.7" "cij=0.8" "cij=0.9" "cij=1"],
       legend=:outerright, 
       ylims=(0,3*10^4),
       xtickfont=font(12), 
      ytickfont=font(12), 
      guidefont=font(12), 
      legendfont=font(12))
plot!(0.01*span_11,Nai_h_c[1,:], color=:black, style = :dash, label="H=0.01")
plot!(0.1*span_11,Nai_h_c[1,:], color=:black, style = :dash, label="H=0.1")
plot!(1*span_11,Nai_h_c[1,:], color=:black, style = :dash, label="H=1")

xlabel!("Exploitation rate (H)")
ylabel!("Adult size (mm)")
png("H_decay_for_cij_gradient_S")




function SLC!(du, u, p, t)
  Na, Sa = u
  t_0, k, r, K, H , d, Smax, gamma, cij = p
   
   
   
   #reproductive cycle
   phi(t) = 2*pi*(t - t_0)/(365)
   periodX(t) = 1/2*(1+tanh(2*sin(phi(t)) - k))


   # reproductive capacity
   avg_size = du[2]
   Smat = 1.34 * (avg_size) - 28.06
   R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax - Smat)), 0.0), 1.0)

  du[1] = dNa = periodX(t)* r * R_ * Na * ((K - Na)/ K) - d * Na - (1 - periodX(t)) * H * Na - cij * Na
  du[2] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - periodX(t))))
end

avg_oocytes = mean([92098, 804183]) # This is the actual mean. mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r_ = reggs*0.998611*0.971057*0.4820525*0.00629
# natural death rates per life stage.
d = 0.000322/365.14  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = 0.00014898749263737575
t0_ = 365.14*0.42
k_= 0.42
Naj_# = 500
K_ = 64000.00    # Carrying capacity
Smax_ = 53.0             # Maximum size for adults3
#         t_0, k,  r,  K,  H ,  d, Smax,  gamma



H1_ = 0

cij_= 1
p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate,cij_,Naj_] 
n=10  #Number of years in the simulation
tspan = (0,365.14*n)
U0_ = [10^4, 49.25]
prob_ = ODEProblem(SLC!, U0_, tspan, p_span)
solve_= solve(prob_, Tsit5())

t_l = length(zeros(Float64,size(1:length(solve_.u))))
n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))
vars = zeros(t_l,n_v)
time2 = zeros(t_l,1)

for j in 1:length(solve_.u)
  for i in 1:2
    vars[j,i] = solve_.u[j][i]
    time2[j] = solve_.t[j]
   end 
end

Na1c0 = vars[:,1]



plot!(Na1c0, label="cij=1",legend=:outerright)

xlims!(0,100)
ylims!(0,7*10^4)
xlabel!("Time (days)", font=12)
ylabel!("Abundance (nº individuals)", font=12)

png("SLC_NAt_H0_cij_3")

#C:\\Users\\miste\\AppData\\Local\\Programs\\Julia-1.9.3\\bin\\julia.exe

SLC_SAt = plot(time2, Sa2, label="Sa SLC")
plot!(time, Sa, label="Sa CLC",xlims=(0,885))
xlabel!("time (days)")
ylabel!("Size (mm)")
png("SLC_SAt")



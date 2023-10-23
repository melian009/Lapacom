
#Packages

using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics


# Simple life cycle on one site.


#Check analytical_SLC_without_S.jl 

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

#=

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
R_ = 500.00
avg_oocytes = 385_613 # This is the actual mean. mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r_= reggs
K_ = 2000    # Carrying capacity


d1_ = 0.590
#H1_ = 0.639
#c12_ = 0.5 #competition term species 2 on 1
Smax_ = 56.0             # Maximum size for adults
Naj = 2000





n=2    #Number of years in the simulation
tspan = length(zeros(Float64,size(1:365.14*n)))
hspan = length(zeros(Float64, size(0:0.1:1)))
span = ones(Float64,size(1:365.14*n))
Kspan = ones(Float64,size(1:365.14*n))*K_ 
R_span = ones(Float64,size(1:365.14*n))*R_
Sm_span = ones(Float64,size(1:365.14*n))*Smax_/2   # Linea de Smax/2 


H_r = range(0, 1, length=hspan)
H_span = ones(Float64,hspan)
cij_span = ones(Float64,hspan)


for i in 1:length(H_span)
    H_span[i] = H_r[i]
    cij_span[i] = H_r[i]
end

H_span
cij_span

Nai = zeros(tspan,hspan)
Sai = zeros(tspan,hspan)
period = zeros(tspan)



for i in 1:length(cij_span)
  cij_ = cij_span[5]
  H_i = H_span[i]
  c=1 
  Nai1 = zeros(Float64,size(1:365.14*n))
  Sai1 = zeros(Float64,size(1:365.14*n))
  periodX = zeros(tspan)

  
  for t_ in 1:(365.14*n)
    t_0 = (365.14*0.42)
    k=0.1
      phi(t_) = 2*pi*(t_ - t_0)/(0.42 * 365)
      periodX[c] = 0.5*(1 - sin(k-phi(t_)))
      
      #Non trivial solution expresions
        Nai1[c] = -(H_i * (1 - periodX[c]) + cij_ * Naj - periodX[c] * r_ * R_ + d1_)* K_ / (periodX[c] * r_ * R_ * 2) 
        Sai1[c] = 1/2*(Smax_ * (1 - H_i * (1 - periodX[c])))
    c=c+1
    end
    Nai[:,i] = Nai1
    Sai[:,i] = Sai1

end

plot(Nai, label=["cij=0" "cij=0.1" "cij=0.2" "cij=0.3" "cij=0.4" "cij=0.5" "cij=0.6" "cij=0.7" "cij=0.8" "cij=0.9" "cij=1"], ylim=(0,K_/2*1.2), legend=:right)
plot!(Kspan/2,label=false, color=:red)
annotate!(120, K_/2*1.03, text("K/2", :red, :center, 7))
ylabel!("Abuncance population")
xlabel!("Time (days)")
plot!(periodX*800, label="X(t)",color=:blue, style = :dash)
plot!(span*800, c=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, 830, text("X(t)=1", :blue, :center, 8))
annotate!(50, 30, text("X(t)=0", :blue, :center, 8))




plot(Sai, label=["H=0" "H=0.1" "H=0.2" "H=0.3" "H=0.4" "H=0.5" "H=0.6" "H=0.7" "H=0.8" "H=0.9" "H=1"], legend=:right)
plot!(Sm_span,color=:red,label = false)
annotate!(-150, Smax_/2, text("Smax/2", :red, :center, 8))

plot!(periodX*10, label="X(t)",color=:blue, style = :dash)
plot!(span*10, color=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, 11, text("X(t)=1", :blue, :center, 8))
annotate!(50, 1, text("X(t)=0", :blue, :center, 8))
xlabel!("Time (days)")
ylabel!("Adult Size (mm)")


#Na and Sa vs H

Sa_H = ones(Float64,hspan)
Na_H = ones(Float64,hspan)

for j in 1:hspan
Sa_H[j] = minimum(Sai[:,j])
Na_H[j] = mean(Sai[:,j])
end
Sa_H[1]

plot(H_span,Sa_H)
plot(H_span,Na_H)










# Non-trivial solution for Na and Sa on a Single Site.
R_ = 50000.00
avg_oocytes = 385_613 # This is the actual mean. mean([92098, 804183])
reggs = avg_oocytes / (365 * 0.42) # conversion rate of adults to eggs.
r_= reggs
K_ = 10000    # Carrying capacity


d1_ = 0.590
#H1_ = 0.639
#c12_ = 0.5 #competition term species 2 on 1
Smax_ = 56.0             # Maximum size for adults
Naj = 2500





n=2    #Number of years in the simulation
tspan = length(zeros(Float64,size(1:365.14*n)))
hspan = length(zeros(Float64, size(0:0.1:1)))
span = ones(Float64,size(1:365.14*n))
Kspan = ones(Float64,size(1:365.14*n))*K_ 
R_span = ones(Float64,size(1:365.14*n))*R_
Sm_span = ones(Float64,size(1:365.14*n))*Smax_/2   # Linea de Smax/2 


H_r = range(0, 1, length=hspan)
H_span = ones(Float64,hspan)
cij_span = ones(Float64,hspan)

for i in 1:length(H_span)
  H_span[i] = H_r[i]
  cij_span[i] = H_r[i]
end



#cij gradient. 

for i in 1:length(H_span)
  cij_ = cij_span[i]
  H_i = H_span[5]
  c=1 
  Nai1 = zeros(Float64,size(1:365.14*n))
  Sai1 = zeros(Float64,size(1:365.14*n))
  periodX = zeros(tspan)

  
  for t_ in 1:(365.14*n)
    t_0 = (365.14*0.42)
    k=0.1
      phi(t_) = 2*pi*(t_ - t_0)/(0.42 * 365)
      periodX[c] = 0.5*(1 - sin(k-phi(t_)))
      
      #Non trivial solution expresions
        Nai1[c] = -(H_i * (1 - periodX[c]) + cij_ * Naj - periodX[c] * r_ * R_ + d1_)* K_ / (periodX[c] * r_ * R_ * 2) 
        Sai1[c] = 1/2*(Smax_ * (1 - H_i * (1 - periodX[c])))
    c=c+1
    end
    Nai[:,i] = Nai1
    Sai[:,i] = Sai1

end



plot(Nai, label=["H=0" "H=0.1" "H=0.2" "H=0.3" "H=0.4" "H=0.5" "H=0.6" "H=0.7" "H=0.8" "H=0.9" "H=1"], ylim=(0,K_/2*1.2), legend=:right)
plot!(Kspan/2,label=false, color=:red)
annotate!(120, K_/2*1.03, text("K/2", :red, :center, 7))
ylabel!("Abuncance population")
xlabel!("Time (days)")
plot!(periodX*800, label="X(t)",color=:blue, style = :dash)
plot!(span*800, c=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, 830, text("X(t)=1", :blue, :center, 8))
annotate!(50, 30, text("X(t)=0", :blue, :center, 8))




plot(Sai, label=["H=0" "H=0.1" "H=0.2" "H=0.3" "H=0.4" "H=0.5" "H=0.6" "H=0.7" "H=0.8" "H=0.9" "H=1"], legend=:right)
plot!(Sm_span,color=:red,label = false)
annotate!(-150, Smax_/2, text("Smax/2", :red, :center, 8))

plot!(periodX*10, label="X(t)",color=:blue, style = :dash)
plot!(span*10, color=:blue, style = :dash, label=false)
plot!(span*0, color=:blue, style = :dash, label=false)
annotate!(50, 11, text("X(t)=1", :blue, :center, 8))
annotate!(50, 1, text("X(t)=0", :blue, :center, 8))
xlabel!("Time (days)")
ylabel!("Adult Size (mm)")

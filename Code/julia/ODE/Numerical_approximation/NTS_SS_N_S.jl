
#Packages

using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots


# Simple life cycle on one site.


#Check analytical_SLC_without_S.jl 

#Exp_lim = 1                 # Exploitation max limit 
#m=0.05                      # Interval of exploitation values 
#for H = 0:m:Exp_lim         #exploitation gradient
#end
function my_ode!(du, u, t, p)
    x1, y1,= u
    re1, K, H1, R1, d1, c12, x2, Smax, g1 = p
    
    function X(t)
      if (t % 365) / 365 >= 0.42
        return 1.0 # Reproductive Cycle
      else
        return 0.0 # Exploitation Cycle
      end
     end
  
    du[1] = (X(t) * re1 * R1 * x1)*((K - x1)/K) - (d1 * x1) - H1*(1 - X(t))*x1 - c12*x1*x2
    du[2] = (g1 * y1) * (1 - (y1/(Smax * (1 - H1 * (1-X(t))))))
    end

#=
function X_(t,t_0,k)
  phi(t) = 2*pi*(t-t_0)/0.42
  X(t) = 0.5*(1*-tanh(k-sen(phi(t))))
  return X(t)
end
=#


u0 = [1000.0, 25.0]  # Initial conditions
tspan = (1.0, 365*2)  # Time span for the simulation (from t=0 to t=1000)
p_1=[0.32, 640000, 0.639, 1, 0.55, 0.5, 10000, 56, 0.6]

prob = ODEProblem(my_ode!, u0, tspan, p_1)
sol = solve(prob, Tsit5()) # "Error: BoundsError: attempt to access Float64 at index [2]"

plot(sol, xlabel="Time", ylabel="State Variables", label=["x1" "x2" "y1" "y2"])



# Non-trivial solution for Na and Sa on a Single Site.

for t_ in 1:365*5 
   function X(t)
    if (t % 365) / 365 >= 0.42
      return 1.0 # Reproductive Cycle
    else
      return 0.0 # Exploitation Cycle
    end
   end

   #Non trivial solution expresions.

   if X(t) == 1
    du[1] = Nai1 = 1/(X*r*R*2*x1)*(K*H*(1-X(t))+c12*x2-X*r*R*K+K*d) #} When X = 1 Reproductive cycle ON.
    du[2] = Sai1 = 1/2*(Smax*(1-H*(1-X)))
   else 
    du[1] = Nai2 = Na0                                             #} When X = 0 Reproductive cycle OFF Nai0 is the last value of abundances from their previous reproductive cycle.
    du[2] = Sai2 = 1/2*(Smax*(1-H*(1-X)))
   end

end


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

R_ = 1.00
r_ = 0.36 #0.32
K_ = 640000.00        # Carrying capacity
d1_ = 0.590
H1_ = 0.639
c12_ = 0.05 #competition term species 2 on 1
Smax_ = 56.0             # Maximum size for adults
Na_0 = 100 # Minimum size of population (lower limit)
p_2 = [R_,r_,K_,d1_,H1_,c12_,Smax_, Na_0]


# Define the ODE problem
prob = ODEProblem(ode_system_solutions!, u0, tspan, p_2)

# Solve the ODE numerically
sol = solve(prob) 
````ERROR: MethodError: no method matching -(::Int64, ::var"#X#8") ```

# Plot the solution
plot(sol, xlabel="Time", ylabel="Na (nº individuals)", label=["Na_i*"])
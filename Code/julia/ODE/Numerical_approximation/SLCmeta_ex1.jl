using DifferentialEquations
using Plots


function ode_system!(du, u, p, t)
    #x1, x2, y1, y2 = u
    #K, r, d = p
    #r = -0.5
    K = 100
    #d = 0.01
   #R1 = 1
   #R2 = 1
   X =  1
   #K = 640000          # Carrying capacity 
   H1 = 0.639
   H2 = 0.57 # Exploitation rate (H)
   #g1 = 0.32
   #g2 = 0.36 # Adult growth rate
    #rates2 = [0.02,0.01]
    #gEA = 0.006          # Instant conversion between stages.
   d1 = 0.59#0.55
   d2 = 0.5    # Natural mortality rate for adults
   #Smax = 56              # Maximum size for adults
   re1 = 0.36 #0.32
   re2 = 0.37 # Adult growth rate
   c12 = 0.05#competition term species 2 on 1
   c21 = 0.025#competition term species 1 on 2
   g1 = 0.36
   g2 = 0.32
   Smax1 = 56
   Smax2 = 50
    
    du[1] = (X * re1 * u[1])*(K - u[1]/K) - (d1 * u[1]) - H1*(1 - X)*u[1] - c12*u[1]*u[2]
    du[2] = (X * re2 * u[2])*(K - u[2]/K) - (d2 * u[2]) - H2*(1 - X)*u[2] - c21*u[1]*u[2]
    du[3] = (g1 * u[3]) * (1 - (u[3]/(Smax1 * (1 - H1 * (1-X)))))
    du[4] = (g2 * u[4]) * (1 - (u[4]/(Smax2 * (1 - H2 * (1-X)))))
    
        
    #du[1] = (r * u[1])*(K - u[1]/K) - d * u[1]
    #du[1] = r * u[1] + u[1]^2 - u[1] * u[2]
    #du[2] = u[1] - u[1] * u[2]
end


u0 = [1.0, 1.0,25.0,25.0]  # Initial conditions for [x, y]
tspan = (0.0, 100.0)  # Time span for the simulation (from t=0 to t=10)


# Define the ODE problem
prob = ODEProblem(ode_system!, u0, tspan)

# Solve the ODE numerically
sol = solve(prob)

# Plot the solution
plot(sol, xlabel="Time", ylabel="State Variables", label=["u[1]" "u[2]" "u[3]" "u[4]"])


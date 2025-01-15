#Packages


using Pkg
Pkg.activate(".")
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics

#=
SLC general model
dNa/dt = (R * re * Na)*((K - Na)/K) - (d * Na) - H*(1 - X(t)) * Na - c12 * Na * x2
dSa/dt = (g * Sa) * (1 - (Sa/(Smax(1 - * H*(1 - X)))))

For X=1 Species reproduces and theres no expoitation on them.
  dNa/dt = (R * re * Na)*((K - Na)/K) - (d * Na) - c12 * Na * x2
  dSa/dt = (g * Sa) * (1 - (Sa/Smax))

Non trivial solution for the scenario of X=1
    
  Nai =  K*(R*re-d-cij*Naj)/(2*R*re)
  Sai =  Snax/2

For X=0 Species reproduces and theres expoitation on them.
  dNa/dt = (R * re * Nai)*((K - Nai)/K) - (d * Na) - H * Na - c12 * Nai * Naj
  dSa/dt = (g * Sa) * (1 - (Sa/Smax))

Non trivial solutuiin for the scenario of X=0
    Nai =  K*(R*re-d-cij*Naj)/(2*R*re)
    Sai = (Smax*(1-H))/2
=#



function SLC!(du, u, p, t)
  Na, Sa = u
  t_0, k, r, K, H , d, Smax, gamma, cij = p
   
   
   
   ## UPDATE for two species 
   
   
   #reproductive cycle
   #phi(t) = 2*pi*(t - t_0)/(365)
   #periodX(t) = 1/2*(1+tanh(2*sin(phi(t)) - k))
   
   function periodX(t)
    if (t % t_0 / t_0) >= k
      return 1.0 # Reproductive Cycle
    else
      return 0.0 # Exploitation Cycle
    end
  end

   # reproductive capacity
   avg_size = du[2]
   Smat = 1.34 * (avg_size) - 28.06
   R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax - Smat)), 0.0), 1.0)

  du[1] = dNa = r * R_ * Na * ((K - Na)/ K) - d * Na - (1 - periodX(t)) * H * Na - cij * Na
  du[2] = dSa = gamma * Sa * (1 - Sa / (Smax - Smax * H * (1 - periodX(t))))
end



## LIST OF PARAMETERS
avg_oocytes = mean([77404, 385613]) # This is the actual mean.
reggs = avg_oocytes / (365 * 0.42) #aplication of the reproduction period stablish by law. (The time banned for extraction or exploitation for the species)
r_ = reggs*0.998611*0.971057*0.4820525*0.00629 # conversion rate of adults to eggs.
# natural death rates per life stage.
d = mean([0.55,0.59])/365.14  # see estimate_mortality_rates.jl for how these values were estimated.
size_growth_rate = mean([0.32,0.36]) #0.00014898749263737575
t0_ = 365.14*0.42
k_= 0.42
Naj_ = 2500
K_ = 64000.00    # Carrying capacity
Smax_ = 53.0             # Maximum size for adults3
#         t_0, k,  r,  K,  H ,  d, Smax,  gamma


n=10    #Number of years in the simulation
t_span = length(zeros(Float64,size(1:365.14*n)))
h_span = length(zeros(Float64, size(0:0.01:1)))
span = ones(Float64,size(1:365.14*n))
#Kspan = ones(Float64,size(1:365.14*n))*K_ 
#Sm_span = ones(Float64,size(1:365.14*n))*Smax_/2   # Linea de Smax/2 


Nai_cij = zeros(t_span,h_span)
Sai_cij = zeros(t_span,h_span)
periodX = zeros(t_span)


H_r = range(0, 1, length=h_span)
H_span = ones(Float64,h_span)
cij_span = ones(Float64,h_span)
#N_H_c_span = zeros(11,11)
#S_H_c_span = zeros(11,11)

for i in 1:length(H_span)
  H_span[i] = H_r[i]
  cij_span[i] = H_r[i]
end

Hs=length(zeros(Float64,size(1:length(H_span))))
cijs = length(zeros(Float64,size(1:length(cij_span))))

#condiciones de estabilidad
H1_ = H_span[1]
cij_= cij_span[1]

p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate,cij_,Naj_]  
n=10  #Number of years in the simulation
tspan = (0,365.14*n)
U0_ = [10^4, 49.25]
prob_ = ODEProblem(SLC!, U0_, tspan, p_span)
solve_= solve(prob_, Tsit5())

t_l = length(zeros(Float64,size(1:length(solve_.u))))
n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))
vars = zeros(t_l,n_v)
time_inicial = zeros(t_l)

for j in 1:length(solve_.u)
  for i in 1:2
    vars[j,i] = solve_.u[j][i]
    time_inicial[j] = solve_.t[j]
  end 
end

#Generacion de metrices cúbicas de almacenamiento 
resultados_simulaciones = zeros(length(time_inicial), length(H_span), length(cij_span)) 
tiempos_totales = zeros(length(time_inicial), length(H_span), length(cij_span))
tiempos_maximos = zeros(length(H_span), length(cij_span))


for i in 1:length(H_span)
  H1_ = H_span[i]
  for j in 1:length(cij_span)

   cij_= cij_span[j]

    p_span = [t0_, k_, r_, K_, H1_, d, Smax_, size_growth_rate,cij_,Naj_]  
    n=10  #Number of years in the simulation
    tspan = (0,365.14*n)
    U0_ = [10^4, 49.25]
    prob_ = ODEProblem(SLC!, U0_, tspan, p_span)
    solve_= solve(prob_, Tsit5())

    t_l = length(zeros(Float64,size(1:length(solve_.u))))
    n_v = length(zeros(Float64,size(1:length(solve_.u[1]))))
    vars = zeros(t_l,n_v)
    time2 = zeros(t_l)
 
      for a in 1:length(solve_.u)
        for b in 1:2
          vars[a,b] = solve_.u[a][b]
          time2[a] = solve_.t[a]
        end 
      end

    if i == 1 && j == 1
      longitud_simulacion = length(time2)
    # Guardar los resultados de la simulación en la matriz cúbica
      resultados_simulaciones[:,i,j] = vars[:,1] # u = vars[:,1] = Abundancias // vars[:,2] = Tallas
      tiempos_totales[:,i,j] = time2
    else 
      lon_0 = length(resultados_simulaciones[:,1,1])
      lon_c = length(vars[:,1])
      elementos_faltantes = lon_0 - lon_c
      if elementos_faltantes > 0 # Concadena el vector corto con zeros para tener la misma longitud del vector de las condiciones iniciales, el mas largo.
        vector_corto = vcat(vars[:,1],zeros(elementos_faltantes))
        vector_corto_t = vcat(time2,ones(elementos_faltantes)*maximum(time2))
      else
        vector_corto = vars[:,1][1:lon_0] # Selecciona los valores del vector mas largo hasta la longitud de la simulación de las condiciones estandar (H=0,c=0)
        vector_corto_t = time2[1:lon_0]
      end
    resultados_simulaciones[:,i,j] = vector_corto
    tiempos_totales[:,i,j] = vector_corto_t
    tiempos_maximos[i,j] = maximum(vector_corto_t) #ultima iteración (nº máximo de días) por simulación
    end
  end
end

resultados_simulaciones
tiempos_totales
tiempos_maximos

#¿En cuáles días se explota y en cuales no (X = 0 y X = 1)? Plot Dynamics

#Initial condition or stable dynamic (H = 0, c = 0)
plot(tiempos_totales[:,1,1],resultados_simulaciones[:,1,1], label="Stable specie")

#to explore exploitation variability without competition (cij = 0)
for i in 1:10:length(H_span)  
  c=1  #c=0.0
  h=i
    if i < length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:blue, label=false)
    else i == length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:blue, label=vcat("cij =",H_span[c]))
    end
end
ylims!(0,6.4*10^4)

# to explore exploitation variability with competition (cij =/= 0)
for i in 1:10:length(H_span)
  c=11 #c_patella_aspera=0.16
  h=i
  if i < length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:green, label=false)
  else i == length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:green, label=vcat("cij =",H_span[c]))
  end
end
ylims!(0,6.4*10^4)

#for patella aspera
for i in 1:10:length(H_span)
  c=17 #c_patella_aspera=0.16
  h=i
  if i < length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:green,style=:dash, label=false)
  else i == length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:green,style=:dash, label=vcat("cij =", H_span[c])) 
  end
end
ylims!(0,6.4*10^4)


for i in 1:10:length(H_span)
  c=51  #cij=0.5
  h=i
  if i < length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:brown, label=false)
  else i == length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:brown, label=vcat("cij =",H_span[c]))
  end
end
ylims!(0,6.4*10^4)

for i in 1:10:length(H_span)
  c=85 #cij=0.84 Patella ordinaria
  h=i
  if i < length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:red, label=false)
  else i == length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:red,style=:dash, label=vcat("cij =",H_span[c]))
  end
end

for i in 1:10:length(H_span)
  c=91 #cij=0.9
  h=i
  if i < length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:red, label=false)
  else i == length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:red, label=vcat("cij =",H_span[c]))
  end
end

for i in 1:10:length(H_span)
  c=101  #cij=1
  h=i
  if i < length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:black, label=false)
  else i == length(H_span)
    plot!(tiempos_totales[:,h,c],resultados_simulaciones[:,h,c],color=:black, label=vcat("cij =",H_span[c]))
  end
end

plot!(tiempos_totales[:,101,1],resultados_simulaciones[:,101,1],color=:black, label=false, legend=:outerright)
xlims!(650,800)
ylims!(2*10^4,6.4*10^4)

#ylims!(4*10^4,7*10^4)
xlims!(152,220) 
#Día es estabilidad, punto medio entre inicio y fin el 
#peródo de explotación X=1 entre los días 612 y 682

#LOCALIZACIÓN DE VALORES PARA UN DÍA ESPECÍFICO
#Matriz de alamacenamiento

time_loc = zeros(Int64,101,101)   #Posición del día a seleccionar dentro de la matrix de tiempos totales-
abun_loc = zeros(Float64,101,101) #Abundancia estimada dentro de la matrix de resultados simulaciones en el día a seleccionar.

dia_estable = 152+(220-152)/2  #Día donde no se aplica la restricción de la restricción de explotación (X=0)

for j in 1:length(H_span)
  for k in 1:length(cij_span)
    for i in 1:length(tiempos_totales[:,1,1])
      if tiempos_totales[i,j,k] > (646) && tiempos_totales[i,j,k] < (648) 
        time_loc[j,k] = i
        abun_loc[j,k] = resultados_simulaciones[i,j,k]
      end
    end
  end
end
time_loc
abun_loc


plot(abun_loc[:,1])
#pendientes!(tiempos_totales, resultados_simulaciones,intervalo,days,comp_Pos)

c=1
days= time_loc[1,1]
X0 = hcat(abun_loc[:,c])
plot(H_span,X0, label=vcat("cij=",cij_span[1]), legend=:outerright)
  

for j in vcat(1,11,17,51,84,91,101) 
  X = abun_loc[:,j]
  intervalo=0.01
  h_span = length(zeros(Float64, size(0:intervalo:1)))
  H_r = range(0, 1, length=h_span)
  H_span = ones(Float64,h_span)

  for i in 1:length(H_span)
    H_span[i] = H_r[i]
  end

  h_N=zeros(length(1:10:length(H_span)))
  h_v=zeros(length(1:10:length(H_span)))
  c=1
  
  for i in 1:10:length(H_span)
    h_v[c] = H_span[i]
    h_N[c] = X[i]
    c=c+1
  end

  h_N_p = vcat(h_N[1],h_N[2],h_N[3],h_N[4],h_N[5],h_N[6],h_N[10],h_N[11])  
  h_v_p = vcat(h_v[1],h_v[2],h_v[3],h_v[4],h_v[5],h_v[6],h_v[10],h_v[11])
  
  # Generación de gráficos coloreados en función del rango de competencia:
  if j == 1 #Non competence
     plot(h_v_p,h_N_p, label=vcat("cij=",cij_span[j]),color=:blue, style=:solid, legend=:bottomleft)#,background=false)
    else
    if j > 1 && j < 12 #10% of competence
     plot!(h_v_p,h_N_p, label=vcat("cij=",cij_span[j]),color=:green, style=:solid)
     else  
     if j > 12 && j < 22 #Patella aspera range of competence
       plot!(h_v_p,h_N_p, label=vcat("c(Patella aspera)=",cij_span[j]),color=:green, style=:dot)
       else
        if j > 22 && j < 52 #50% of competence
          plot!(h_v_p,h_N_p, label=vcat("cij=",cij_span[j]),color=:red, style=:solid)
          else 
          if j > 52 && j < 85 #Patella ordinaria range of copetence
            plot!(h_v_p,h_N_p, label=vcat("c(Patella ordinaria)=",cij_span[j]),color=:brown, style=:dot)
            else 
            if j > 85 && j < 92 #90% of competence
              plot!(h_v_p,h_N_p, label=vcat("cij=",cij_span[j]),color=:brown, style=:solid)
              else j > 93 && j < 102 #100% of competence
                plot!(h_v_p,h_N_p, label=vcat("cij=",cij_span[j]),color=:black, style=:solid)
            end
          end
        end
      end
    end
  end
end
xlabel!("Exploitation rate (H)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("Exploitation vs Abuncance")
savefig("FIG_4a_legend.png")



#=
plot(H_span,X0[:,2],label="cij=0.0", color=:blue)
plot!(H_span,X1[:,2],label="cij=0.1", color=:green)
plot!(H_span,pa[:,2],label="c(pa,po)=0.16", color=:green, style = :dash)
plot!(H_span,X2[:,2],label="cij=0.5", color=:red)
plot!(H_span,po[:,2],label="c(po,pa)=0.83", color=:brown, style = :dash)
plot!(H_span,X3[:,2],label="cij=0.9", color=:brown)
plot!(H_span,X4[:,2],label="cij=1.0", color=:black,legend=:bottomleft,
      background=nothing)
xlabel!("Exploitation rate (H)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("X=0")
savefig("SLC_NA_H_cij_X0.png")


plot!(H_span,X6[:,2],label="cij=0.1", color=:green,style=:solid)
plot!(H_span,X7[:,2],label="cij=0.5", color=:red, style=:solid)
plot!(H_span,X8[:,2],label="cij=0.9", color=:brown, style=:solid)
plot!(H_span,X9[:,2],label="cij=1.0", color=:black, style=:solid)
xlabel!("Exploitation rate (H)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("X=1")
savefig("SLC_NA_H_cij_X1.png")
=#


#H=0.0
plot(tiempos_totales[:, 1,1],resultados_simulaciones[:, 1,1], label="cij=0.00", color=:blue, style=:solid)
plot!(tiempos_totales[:,1,11],resultados_simulaciones[:, 1,11],label="cij=0.10",  color=:blue, style=:dash)
plot!(tiempos_totales[:,1,51],resultados_simulaciones[:, 1,51],label="cij=0.50",  color=:blue, style=:dashdot)
plot!(tiempos_totales[:,1,91],resultados_simulaciones[:, 1,91],label="cij=0.90", color=:blue, style=:dashdotdot,
      background=nothing)

xlims!(0,365.14*8)
ylims!(5*10^4,7*10^4)
xlabel!("Time (days)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("H=0.0")
savefig("Figure_4c.png")

#H=0.5
plot(tiempos_totales[:, 51,1],resultados_simulaciones[:, 51,1], label="cij=0.00", color=:green, style=:solid)
plot!(tiempos_totales[:, 51,11],resultados_simulaciones[:, 51,11],label="cij=0.10",  color=:green, style=:dash)
plot!(tiempos_totales[:, 51,51],resultados_simulaciones[:, 51,51],label="cij=0.50",  color=:green, style=:dashdot)
plot!(tiempos_totales[:, 51,91],resultados_simulaciones[:, 51,91],label="cij=0.90", color=:green, style=:dashdotdot,
      background=nothing)

xlims!(0,365.14*8)
ylims!(5*10^4,7*10^4)
xlabel!("Time (days)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("H=0.5")
savefig("Figure_4d.png")

#H=0.9
plot(tiempos_totales[:, 91,1],resultados_simulaciones[:, 91,1], label="cij=0.00", color=:red, style=:solid)
plot!(tiempos_totales[:, 91,11],resultados_simulaciones[:, 91,11],label="cij=0.10",  color=:red, style=:dash)
plot!(tiempos_totales[:, 91,51],resultados_simulaciones[:, 91,51],label="cij=0.50",  color=:red, style=:dashdot)
plot!(tiempos_totales[:, 91,91],resultados_simulaciones[:, 91,91],label="cij=0.90", color=:red, style=:dashdotdot,
      background=nothing)

xlims!(0,365.14*8)
ylims!(5*10^4,7*10^4)
xlabel!("Time (days)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
title!("H=0.9")
savefig("Figure_4e.png")





#Simulations for Patella ordinaria and Patella asoera 
oocytes_po = 385613 # Average: Patella ordinaria (nº of Eggs)
oocytes_pa = 77404  # Average: Patella aspera (nº of Eggs)
c_po = oocytes_po/(oocytes_pa + oocytes_po)
c_pa = oocytes_pa/(oocytes_po + oocytes_pa)
oocytes = [oocytes_po,oocytes_pa]    # Patella ordinaria, Patella aspera
reggs = oocytes / (365 * 0.42)       # Population growth rate    
re = reggs*0.998611*0.971057*0.4820525*0.00629     # conversion of poplation growth rate from eggs to adults.
Kt = 64000          # Carrying capacity
H_ = [0.639,0.57] # Exploitation rate (H)
da_ = [0.55,0.59]    # Natural mortality rate for adults
Sm = 56              # Maximum size for adults
gammas = [0.32,0.36] # Adult growth rate
i = [1,2]            # Species: "Patella ordinaria" (i=1); "Patella aspera" (i=2)
Naj_po = reggs[1]
Naj_pa = reggs[2]
k_ = 0.42
t0_ = 365.25

# Patella aspera simulation
p_span_po = [t0_, k_, reggs[1], K_, H_[1], da_[1], Smax_, gammas[1],c_pa,Naj_pa]  
n=10  #Number of years in the simulation
tspan = (0,365.14*n)
U0_ = [10^3, 49.25]
prob_ = ODEProblem(SLC!, U0_, tspan, p_span_po)
solve_= solve(prob_, Tsit5())

t_l_po = length(zeros(Float64,size(1:length(solve_.u))))
n_v_po = length(zeros(Float64,size(1:length(solve_.u[1]))))
vars_po = zeros(t_l_po,n_v_po)
time2_po = zeros(t_l_po,1)

for j in 1:length(solve_.u)
    for i in 1:2
      vars_po[j,i] = solve_.u[j][i]
      time2_po[j] = solve_.t[j]
    end 
end
vars_po
time2_po

# Patella aspera simulation
p_span_pa = [t0_, k_, reggs[2], K_, H_[2], da_[2], Smax_, gammas[2],c_po,Naj_po]  
n=10  #Number of years in the simulation
tspan = (0,365.14*n)
U0_ = [10^4, 49.25]
prob_ = ODEProblem(SLC!, U0_, tspan, p_span_pa)
solve_= solve(prob_, Tsit5())

t_l_pa = length(zeros(Float64,size(1:length(solve_.u))))
n_v_pa = length(zeros(Float64,size(1:length(solve_.u[1]))))
vars_pa = zeros(t_l_pa,n_v_pa)
time2_pa = zeros(t_l_pa,1)

for j in 1:length(solve_.u)
    for i in 1:2
      vars_pa[j,i] = solve_.u[j][i]
      time2_pa[j] = solve_.t[j]
    end 
end

#Patella ordinaria outputs
vars_po
time2_po
#Patella aspera outputs
vars_pa
time2_pa


#Igualar longitud de vectores para representar ambas señales en una misma gráfica

#if length(vars_po[:,1]) > length(vars_pa[:,1])
  mat_spp = zeros(length(time2_po),2) 
  mat_spp_t = zeros(length(time2_po),2)
  longitud_simulacion = length(time2_po)
# Guardar los resultados de la simulación en la matriz cúbica
  mat_spp[:,1] = vars_po[:,1] # u = vars[:,1] = Abundancias // vars[:,2] = Tallas
  mat_spp_t[:,1] = time2_po

  lon_0 = length(mat_spp[:,1])
  lon_c = length(vars_pa[:,1])
  elementos_faltantes = lon_0 - lon_c
  vector_corto = vcat(vars_pa[:,1],zeros(elementos_faltantes))
  vector_corto_t = vcat(time2_pa,ones(elementos_faltantes)*maximum(time2_pa))
  mat_spp[:,2] = vector_corto
  mat_spp_t[:,2] = vector_corto_t
#=else
  mat_spp = zeros(length(time2_pa),2) 
  mat_spp_t = zeros(length(time2_pa),2)
  longitud_simulacion = length(time2_pa)
# Guardar los resultados de la simulación en la matriz cúbica
  mat_spp[:,1] = vars_pa[:,1] # u = vars[:,1] = Abundancias // vars[:,2] = Tallas
  mat_spp_t[:,1] = time2_pa

  lon_0 = length(mat_spp[:,1])
  lon_c = length(vars_po[:,1])
  elementos_faltantes = lon_0 - lon_c
  vector_corto = vcat(vars_po[:,1],zeros(elementos_faltantes))
  vector_corto_t = vcat(time2_po,ones(elementos_faltantes)*maximum(time2_po))
  mat_spp[:,2] = vector_corto
  mat_spp_t[:,2] = vector_corto_t
# end
=#

#Plotting species

plot(mat_spp_t[:,1],mat_spp[:,1],label="Patella ordinaria",  color=:blue, style=:solid)
plot!(mat_spp_t[:,2],mat_spp[:,2],label="Patella aspera",  color=:orange, style=:solid, background=nothing)
xlims!(000,2000)
ylims!(3*10^4,6.4*10^4)
xlabel!("Time (days)", font=12)
ylabel!("Abundance (nº individuals)", font=12)


plot(tiempos_totales[:,64,85],resultados_simulaciones[:,64,85],label="Patella aspera",  color=:orange, style=:solid)
plot!(tiempos_totales[:,57,17],resultados_simulaciones[:,57,17],label="Patella ordinaria",  color=:blue, style=:solid, background=nothing)
xlims!(000,2000)
xlabel!("Time (days)", font=12)
ylabel!("Abundance (nº individuals)", font=12)
savefig("Figure_4b.png")



contour(tiempos_totales)







# Non-trivial solution for Na and Sa on a Single Site.
abun_loc #this matrix is filled by the average abundance obtained by a numerical aproximation

#The exercise here is to obtain the same or at least a similar matrix from analitical aproximation.
#As more similar, more presition will have the model.
# Contrast between hexploitation traits on each competition values.
#Parameters and initial conditions

avg_oocytes = [77404, 385613] # This is the actual mean. 
reggs = avg_oocytes / (365.14 * 0.42) # conversion rate of adults to eggs.
r = reggs*0.998611*0.971057*0.4820525*0.00629 # conversion rate of adults to eggs.
r_ = mean(r)

K_ = 64000    # Carrying capacity
d = [0.55,0.59]/365.14 #Valores individuales  (día^{-1})
d_= mean(d) #valor promedio
Smax_ = 56.0             # Maximum size for adults
Naj = mean(r)
Smat = mean([[33.4,37.4],[34.6,37.5]])
Smat_ = mean(Smat)
avg_size_0 = 33.4  #Average size as mean size at first maturity as an initial condition



#Vectors for ploting and simulations
n=10   #Number of years in the simulation
t_span = length(zeros(Float64,size(1:365.14*n))) #Number of years in days for the analytical sumulations
v_span = length(zeros(Float64, size(0:0.01:1)))   #Number of iintervals that our varables to contrast will be considered
span = ones(Float64,size(1:365.14*n))

#Kspan = ones(Float64,size(1:365.14*n))*K_ 
#Sm_span = ones(Float64,size(1:365.14*n))*Smax_/2   # Linea de Smax/2 

#Storage exploitation matrix for a specific competition values cij=constant
Nai_h = zeros(t_span,v_span)
Sai_h = zeros(t_span,v_span)
periodX = zeros(t_span)

#Rangos de las variables a contrastar
H_r = range(0, 1, length=v_span)
t_r = range(0, t_span, length=t_span)

#Theoretical values for exploitation and competiticon 

H_span = ones(Float64,v_span)
cij_span = ones(Float64,v_span)

for i in 1:length(H_span)
   
  H_span[i] = H_r[i]
  cij_span[i] = H_r[i]

end
#Outputs
H_span
cij_span


# Time values for analitical stimation

vector_t = ones(Float64,t_span)

for i in 1:length(vector_t)
 
  vector_t[i] = t_r[i]

end
#Output
vector_t

#Outputs of analitical aproximation
Nai_h = zeros(Float64,length(H_span),length(cij_span))
Sai_h = zeros(Float64,length(H_span),length(cij_span))


# Analitical aproximation for Exploited scenario on the metacommunity dynamic0
for j in 1:length(cij_span) 
  cij_= cij_span[j]

  Nai1=zeros(Float64,length(vector_t))
  Sai1=zeros(Float64,length(vector_t))
    
 for  i in 1:length(H_span)
  H_i = H_span[i]
    avg_size = 49.25
    k=0.42
    t_0 = 365.24 
    for t_ in 1:length(vector_t)
      if (t_ % t_0)/t_0 >= k_
        X=1.0
      else
        X=0.0
      end

    Smat = 1.34 * (avg_size) - 28.06
    R_ = min(max(0.5 * (1.0 + (avg_size - Smat) / (Smax_ - Smat)), 0.0), 1.0)
    #Nom trivial solution for Exploitation scenario (X=0)
      #Abundance
      Nai1[t_] = (K_- ((-H_i -cij_ * Naj - d_+ r_ * R_)*(r_ * R_ * 2)/K_))
      #Adult Size
      Sai1[t_] = (Smax_ * (1 - H_i))/2
      avg_size = Sai1[t_]
    end

    Nai_h[i,j] = mean(Nai1)
    Sai_h[i,j] = mean(Sai1)
 end
end

   
Nai_h

plot(H_span,Nai_h[:,1])

for i in 1:10:length(H_span)  
    if i ==1
    plot(H_span,Nai_h[:,i],color=:blue, label=vcat("cij =",cij_span[i]))  
    else
      if i < length(H_span)
      plot!(H_span,Nai_h[:,i],color=:blue, label=vcat("cij =",cij_span[i]))
      else i == length(H_span)
      plot!(H_span,Nai_h[:,i],color=:blue, label=vcat("cij =",cij_span[i]))
      end
    end
end
xlims!(0.0,1)
ylims!(0,6.4*10^4)

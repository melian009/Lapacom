The graphics in this folder refer to the dynamics represented in the file NTS_SS_N_S.jl. Here's a 
breakdown:
-	HX_NAt.png:
   This figure represents the temporal simulation of the non-trivial solution of the SLC model 
(Eq: 3 from Overleaf). The legend shows the exploitation values for which each line 
(population) is represented. X(t) represents the reproductive cycle of limpets, defined as a 
smoothed periodic function. Except for the exploitation value, the rest of the parameters used 
remain the same for all lines. The maximum population when H=0 is K/2, where K is the 
carrying capacity of the species.
-	Hx_SAt:
   This figure represents the same simulation as HX_NAt.png, but with adult size instead of 
abundance. As exploitation increases when the reproductive cycle ends (X(t) ? 0), the 
difference in exploitation levels increases.
-	cijX_NAt:
   This figure represents the same situation as HX_Nat.png, but instead of varying the 
exploitation coefficient, the model fixes the exploitation at 0.6 and calculates dynamics for 
different cij values. Cij is a parameter measuring the spatial interaction a species has in the 
model, represented as cij xi xj, where i and j are species competing for the same space, and cij 
is the normalized interaction or competition parameter between 0 and 1.
-	cijX_SAt:
   This figure represents the same simulation as cijX_NAt.png, but with adult size instead of 
abundance. It shows that the size of the models depends only on the degree of exploitation 
and not on other variables that may affect abundance.
H_Na and H_Sa:
      These figures represent the exploitation gradient in abundance and adult size. The 
maximum abundance and minimum size of adults are plotted against the different simulated 
exploitation values (from 0 to 1). These values (Na max and Sa min) were chosen because the 
real values of the minimum abundance and maximum size of the represented dynamics are 0 
and Smax/2, respectively.
-	Cij_Na and cij_Sa:
   Same situation as the previous point but against cij instead of H. Abundance shows a 
negative slope as the competition effect increases, but it is constant in size.
-	CLC_Nat:
   Representation of the egg and trochophore stages in the 10-year CLC dynamics on a 
vertical scale between 10^5 and 4*10^7 individuals and a horizontal scale of 1000 days.
-	CLC_Nat2:
   Representation of the egg, trochophore, veliger, juveniles, and adults stages in the CLC 
dynamics. The vertical scale is adapted to visualize all stages below 10^5 individuals, and the 
time scale is the same as for CLC_NAt.
-	CLC:
   This figure is a composition of the two previous ones, visualizing the complete dynamics of 
the CLC for each stage.
-	CLC_SA:
   Representation of the evolution of adult size over time for the CLC simulated with an 
exploitation of H = 0.6.


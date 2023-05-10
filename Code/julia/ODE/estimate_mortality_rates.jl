# I know the mortality rate of the limpet is 0.55 per year but this number is across all life stages. My goal is to find a stage-specific mortality rate to parametarize my model. Another thing we know is that mortality rate decreases sharply as the limpet goes to later life stages.I know the mortality rate of the limpet is 0.55 per year but this number is across all life stages. My goal is to find a stage-specific mortality rate to parametarize my model. Another thing we know is that mortality rate decreases sharply as the limpet goes to later life stages.
# I have estimated the mortality rate of adults based on the fact that they live 7 to 10 years. It's   1 / age. 1/3102 = 0.000322 per day (time step of the model) or 0.1175 per year. 3102 is about 8.5 years.
# I would like to use fit a exponential decay function such that it has mean of 0.55 and mortality rate of adults is 0.1175.

# First, let's define the life stages numerically: 0 for egg, 1 for trocophore, 2 for veliger, 3 for juvenile, and 4 for adult.
# exponential decay function:
# Mortality_rate(stage) = A * exp(-B * (4 - stage))
# Here, "A" is the mortality rate at the adult stage, "B" is the decay constant, and "stage" represents the life stage. This equation ensures that the mortality rate is highest at the egg stage (stage = 0) and decreases as the limpet progresses through its life stages.


using Optim

adult_mortality_rate = 0.1175  # P. ordinaria

# Initial guess for B
B0 = [0.1]

# Define the life stages
stages = [1, 2, 3, 4, 5]  # Enumerating life stages starting from 1 (Egg, Trocophore, Veliger, Juvenile, Adult)

# Define the exponential decay function
exp_decay(stage, A, B) = A .* exp(-B .* (5 .- stage))

# Objective function to minimize
function objective(B)
  mortality_rates = exp_decay.(stages, adult_mortality_rate, B[1])
  mean_mortality_rate = sum(mortality_rates) / length(stages)
  return (mean_mortality_rate - 0.55)^2
end


# Minimize the objective function
result = optimize(objective, B0)

# Get the optimal value of B
B_opt = Optim.minimizer  # = -0.6026855468750002

# Death rates: 
# Stage d
# 1 1.309211754648921
# 2 0.7165836397346279
# 3 0.3922147130988946
# 4 0.2146747045860786
# 5 0.1175


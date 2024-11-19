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
  # Check if any mortality rate is greater than 1
  for rate in mortality_rates
    if rate > 1
      return Inf  # Return a large penalty
    end
  end
  mean_mortality_rate = sum(mortality_rates) / length(stages)
  return (mean_mortality_rate - 0.55)^2
end


# Minimize the objective function
result = optimize(objective, B0)

# Get the optimal value of B
B_opt = Optim.minimizer  # = -0.5353292301297188

# Death rates: 
# Stage d
# 1 0.999999975120952
# 2 0.585476487339919
# 3 0.34278272575599816
# 4 0.20069123118943133
# 5 0.1175

# theoretically, there should only be a single value of B that minimizes the function and meets the constraints of maintaining the mean mortality rate at 0.55. This is because the exponential function is strictly decreasing (if B > 0) or increasing (if B < 0), so it will not have local minima or maxima -- only a global minimum or maximum.
# In our specific problem, we are minimizing the squared difference between the mean mortality rate and the target value (0.55). This objective function forms a parabola in the B-dimension, and parabolas have a unique minimum (in the case of a convex parabola, which is our case here).
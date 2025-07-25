datadir = "../../../Code/julia/ODE/"
nsites = 8
nlifestages = 5

species = ["P. ordinaria", "P. aspera"]
nspecies = length(species)

## initial state of the system
"""
    initialize_species(nspecies, nsites, nlifestages)

Initialize the population and average size for a given number of species, sites, and life stages.

# Arguments
- `nspecies`: Number of species.
- `nsites`: Number of sites.
- `nlifestages`: Number of life stages.

# Returns
- `u0_general_all`: A 3D array where the first dimension represents sites, the second dimension represents life stages plus average size, and the third dimension represents species. Each element in the array represents the population of a life stage (or the average size) at a site for a species.

# Example
```julia
u0_general = initialize_species(2, 10, 5)

"""
function initialize_species(nspecies, nsites, nlifestages; initial_population=[58_000.0, 58_000.0], initial_size=[48.0, 48.0])
  u0_general_all = []
  for species in 1:nspecies
    # Initialize for each species
    u0_general = [[initial_population[species] for j in 1:nlifestages] for i in 1:nsites]
    # Initial average sizes per site. Avg. size is between 45 to 51.
    for i in 1:nsites
      push!(u0_general[i], initial_size[species])
    end
    # u0 cannot be a nested vector.
    u0_general = reshape(reduce(vcat, u0_general), length(u0_general[1]), length(u0_general))'
    push!(u0_general_all, u0_general)
  end
  # Combine all species into a 3D array
  u0_general_all = cat(u0_general_all..., dims=3)
  return u0_general_all
end

u0_general = initialize_species(2, nsites, nlifestages)

## defining model parameters

"""
captures the reproductive cycle and equals to zero when fishing and to one when organisms reproduce. Reproduction occurs between November and March (0.42 of the year).
"""
function reproductive_cycle(t)
  if (t % 365) / 365 >= 0.42
    return 0.0
  else
    return 1.0
  end
end

# conversion rates between stages
# average oocytes per year per adult
avg_oocytes_1 = 385_613 # This is the actual mean. mean([92098, 804183])
reggs_1 = avg_oocytes_1 / (365 * 0.42) # conversion rate of adults to eggs.
r_1 = [reggs_1, 0.998611, 0.971057, 0.4820525, 0.00629]

avg_oocytes_2 = 77_101 # This is the actual mean. mean([92098, 804183])
reggs_2 = avg_oocytes_2 / (365 * 0.42) # conversion rate of adults to eggs.
r_2 = [reggs_2, 0.998611, 0.971057, 0.4820525, 0.00629]

r = cat(r_1, r_2, dims=2)

# natural death rates per life stage.
function return_death_rate()
  d = [0.999 / 365, 0.585 / 365, 0.343 / 365, 0.201 / 365, 0.000322]  # see estimate_mortality_rates.jl for how these values were estimated.
  for species in 2:nspecies
    d = cat(d, d, dims=2)
  end
  return d
end

d = return_death_rate()
size_growth_rate = [0.00014898749263737575, 0.00014898749263737575]

distance_df = CSV.read(joinpath(datadir, "distance_matrix.csv"), DataFrame)
distance_matrix = Float64.(Matrix(distance_df)[:, 1:end-1])
distance_matrix = distance_matrix[1:nsites, 1:nsites]
migraation_probs_df = CSV.read(joinpath(datadir, "migration_probabilites_among_sites.csv"), DataFrame)
# Modify distances (migration probs) by movement prob. which are rough estimates based on the land shape and oceanic currents
mig_probs = Float64.(Matrix(migraation_probs_df)[:, 1:end-1])
mig_probs = mig_probs[1:nsites, 1:nsites]


# exploitation_rates = [0.436795998061044, 0.43767209156532155, 0.4603254055329175, 0.38748435327632225, 0.40337922500828566, 0.5105131417482706, 0.4799123913754184, 0.47959950031955256]  # These values are estimated based on the average size of the limpets at each site. See the text below.
exploitation_rates = [0.9880685477383748, 0.9880947574980841, 0.9891331673866427, 0.9874718910553065, 0.9875900832783625, 0.9900632400512024, 0.9894700088322096, 0.9894645471975891] # These values are estimated based on the average size of the limpets at each site. See the text below. These second exploitation rates used Gompertz growth function instead of logistic.
exploitation_rates = exploitation_rates .- exploitation_rates[4]  # subtract the exploitation rate of Desertas from all sites because it is a fully protected area and no fishing happens there. This is to make sure that the exploitation rate of Desertas is zero and the exploitation rates of other sites are relative to Desertas.

exploitation_rates = cat(exploitation_rates, exploitation_rates, dims=2) # for two species

# function exploitation_rate(size, site)
#   size_max = 56.0
#   min_size_for_fishing = 40
#   exploitation_rates_max = [0.43679599501684097, 0.43767209022758713, 0.46032540377965536, 0.38748435785424035, 0.4033792240083418, 0.5105131408815863, 0.47991239051754064, 0.4795994994034682]  # These values are estimated based on the average size of the limpets at each site.
#   # exploitation_rates_max = [0.980700696804078, 0.9807285958927162, 0.989105147959682, 0.9796381668584575, 0.9875759368524103, 0.9900648704172368, 0.9894851154977173, 0.9894795735668236] # These exploitation rates used Gompertz growth function instead of logistic.
#   exploitation_rates_max = exploitation_rates_max .- exploitation_rates_max[4]  # subtract the exploitation rate of Desertas from all sites because it is a fully protected area and no fishing happens there. This is to make sure that the exploitation rate of Desertas is zero and the exploitation rates of other sites are relative to Desertas.
#   if size <= min_size_for_fishing
#     return 0.0
#   else
#     slope = 1 / (size_max - min_size_for_fishing)
#     intercept = -min_size_for_fishing / (size_max - min_size_for_fishing)
#     E = slope * size + intercept
#     E = min(E, 1.0)
#     return E * exploitation_rates_max[site]
#   end
# end

size_max = [56.0, 56.0]
K = [64_000, 64_000]  # for 6.4 km2 per site.
α = [0.1, 0.1, 0.1]  # dispersion factor for Egg, Trochophore, and Veliger
α = cat(α, α, dims=2) # for two species
# Since we do not have any info about site size, dispersion is only a function of dispersion factor and distance.

competition_coefficient = [0.5, 0.5]  # first element is the effect of competition of species 1 on species 2 and the second element is the effect of competition of species 2 on species 1.

p_general = [r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α, mig_probs, competition_coefficient]

"""
I flatten the parameters and then return them back in the original shape inside `nsites!`. This is because to run Sensitivity analysis (`gsa`), I have to pass flat bounds list.
"""
function flatten_p(p)
  flat_p = Real[]
  for p_list in p
    if typeof(p_list) <: Real
      push!(flat_p, p_list)
    else
      for p in p_list
        push!(flat_p, p)
      end
    end
  end
  return Tuple(flat_p)
end

p_general_flat = flatten_p(p_general)

function restructure_flat_p(flat_params::Tuple, original_shape)
  flat_params = [i for i in flat_params]
  nested_params = []
  start_idx = 1
  for p_list in original_shape
    if typeof(p_list) <: Real
      push!(nested_params, flat_params[start_idx])
      start_idx += 1
    else
      end_idx = start_idx + length(p_list) - 1
      push!(nested_params, reshape(eltype(p_list).(flat_params[start_idx:end_idx]), size(p_list)))
      start_idx = end_idx + 1
    end
  end
  return Tuple(nested_params)
end

"""
  Estimates the `size at first maturity` by constructing a linear line between `maximum size` and average size before and after protection (These values come from empirical data).

  Reproduction capacity depends on size where at maximum size, reproduction is 100% and at size at first maturity it is 50%.

  NOTE The function only implements the calculations for P. ordinaria. Write a new function fo P. aspera.

  ## Data


  ### Size at first maturity

  Size at first Maturity bef/after [P. ordinaria,P. aspera]  = [33.4/37.4, 34.6/37.5]

  ### Average size

  Before (1996-2006): FULL ACCESS.  # Use this for before
  Patella apera = 43.53mm
  Patella ordinaria = 46.26mm

  Only MPA  # Use this for after
  After (2007-2017)
  Patella aspera = 50.61mm
  Patella ordinaria = 49.25mm

  ## The model

  Let's represent size at first maturity by "M" and the average body size by "A". To establish the relationship between the two sizes, we can use a simple linear equation:

  M = kA + b

  where "k" is the proportionality constant and "b" is a constant term.

  We have two sets of values, one for the time before fishing protection (M₁/A₁ = 33.4/46.26) and one for the time after fishing protection (M₂/A₂ = 37.4/49.25). To find the values of "k" and "b", we can set up two equations:

  33.4 = k * 46.26 + b
  37.4 = k * 49.25 + b
  Solving this system of linear equations gives us values of "k" and "b".

  M ≈ 1.34A - 28.06
"""
function calculate_size_at_first_maturity_ordinaria(current_avg_size)
  M = 1.34 * (current_avg_size) - 28.06
end

"""
  ## The model

  Let's represent size at first maturity by "M" and the average body size by "A". To establish the relationship between the two sizes, we can use a simple linear equation:

  M = kA + b

  where "k" is the proportionality constant and "b" is a constant term.

  We have two sets of values, one for the time before fishing protection (M₁/A₁ = 34.6/43.53) and one for the time after fishing protection (M₂/A₂ = 37.5/50.61). To find the values of "k" and "b", we can set up two equations:

  34.6 = k * 43.53 + b
  37.5 = k * 50.61 + b
  Solving this system of linear equations gives us values of "k" and "b".

  M ≈ 0.4096045 A - 16.769916
"""
function calculate_size_at_first_maturity_aspera(current_avg_size)
  M = 0.4096045 * (current_avg_size) - 16.769916
end

calculate_size_at_first_maturity(current_avg_size, species) = species == 1 ? calculate_size_at_first_maturity_ordinaria(current_avg_size) : calculate_size_at_first_maturity_aspera(current_avg_size)

"""
  Return the reproduction capacity (between 0 and 1) given the current average size and size at first maturity and maximum size
"""

"Return the reproduction capacity (between 0 and 1) given the current average size and size at first maturity and maximum size"
reproduction_capacity(Saverage, Smaturity, Smax) = min(max(0.5 * (1.0 + (Saverage - Smaturity) / (Smax - Smaturity)), 0.0), 1.0)

function nsites!(du, u, p, t)
  original_shape = ([2515.414220482714 502.94194390084806; 0.998611 0.998611; 0.971057 0.971057; 0.4820525 0.4820525; 0.00629 0.00629], [0.002736986301369863 0.002736986301369863; 0.0016027397260273972 0.0016027397260273972; 0.0009397260273972603 0.0009397260273972603; 0.0005506849315068494 0.0005506849315068494; 0.000322 0.000322], [0.00014898749263737575, 0.00014898749263737575], [0.0 12.842542485919632 33.20921059959618 73.75112848881275 45.78816486689583 40.95386826936409 23.127098263365838 13.53847480618783; 12.842542485919632 0.0 31.551937191782326 74.4358482251258 50.57366923475108 43.525278937384996 18.68655292046152 19.67083747567315; 33.20921059959618 31.551937191782326 0.0 43.012716906931786 25.586975451710284 15.767470469506062 13.16139561432816 21.70783981938122; 73.75112848881275 74.43584822512581 43.012716906931786 0.0 32.56225153468865 33.114636105203594 56.164462942582226 60.38662955616323; 45.78816486689583 50.573669234751094 25.586975451710284 32.56225153468865 0.0 10.08940978792187 35.940503401029886 32.44536267645653; 40.95386826936408 43.525278937384996 15.767470469506062 33.114636105203594 10.08940978792187 0.0 27.21061174752855 27.454025172124666; 23.127098263365838 18.68655292046152 13.16139561432816 56.164462942582226 35.940503401029886 27.21061174752855 0.0 16.029653650497625; 13.53847480618783 19.67083747567315 21.70783981938122 60.38662955616323 32.44536267645653 27.454025172124666 16.029653650497625 0.0], [0.4940342738691874 0.4940342738691874; 0.49404737874904203 0.49404737874904203; 0.49456658369332135 0.49456658369332135; 0.49373594552765326 0.49373594552765326; 0.49379504163918125 0.49379504163918125; 0.4950316200256012 0.4950316200256012; 0.4947350044161048 0.4947350044161048; 0.49473227359879457 0.49473227359879457], [56.0, 56.0], [64000, 64000], [0.1 0.1; 0.1 0.1; 0.1 0.1], [0.0 0.8 0.0 0.0 0.0 0.0 0.2 0.8; 0.0 0.0 0.2 0.0 0.0 0.2 0.5 0.0; 0.0 0.6 0.0 0.4 0.2 0.4 0.8 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.2 0.8 0.8 0.0 1.0 0.6 0.0; 0.0 0.4 0.8 0.6 0.2 0.0 0.6 0.0; 0.0 0.8 0.4 0.2 0.2 0.4 0.0 0.0; 0.8 0.6 0.0 0.0 0.0 0.0 0.0 0.0], [0.5, 0.5])
  p = restructure_flat_p(p, original_shape)
  nsites = 8
  nlifestages = 5
  reduce_trochophore_dispersals = 0.2 # dispersal rates of Trochophores should be a fraction of Eggs. This factor is by how much it is reduced.

  r, d, size_growth_rate, distance_matrix, exploitation_rates, size_max, K, α, mig_probs = p
  for site in 1:nsites

    for species in 1:nspecies
    
      # Stages 1 Egg
      stage = 1
      prev_stage = 5

      # immigration
      dispersal_probs1 = (exp(-α[stage, species]) ./ distance_matrix[:, site])
      dispersal_probs1[site] = 0.0
      dispersal_probs1 = dispersal_probs1 .* mig_probs[:, site]
      # dispersal_probs1 = dispersal_probs1 ./ sum(dispersal_probs1)

      # emigration
      dispersal_probs2 = (exp(-α[stage, species]) ./ distance_matrix[site, :])
      dispersal_probs2[site] = 0.0
      dispersal_probs2 = dispersal_probs2 .* mig_probs[site, :]
      # dispersal_probs2 = dispersal_probs2 ./ sum(dispersal_probs2)

      Saverage = du[site, 6, species]  # 6 is the index of average size
      Smaturity = calculate_size_at_first_maturity(Saverage, species)

      du[site, stage, species] = (reproductive_cycle(t) * r[stage, species] * u[site, prev_stage, species] * reproduction_capacity(Saverage, Smaturity, size_max[species])) -
                        (r[stage+1, species] * u[site, stage, species]) -
                        (d[stage, species] * u[site, stage, species]) +
                        (sum(dispersal_probs1 .* u[:, stage, species])) -
                        (u[site, stage, species] * sum(dispersal_probs2))
      # What does the last part of the first term do? We use it to reduce the reproductive capacity of the adults based on their size. At max size, 100% of them reproduce. At "size at first maturity" size, only 50% of them reproduce. Make a ramp function that calculates the fraction of the adults that reproduce given the mean size of the individuals. For this, we will also need to have another equation (`calculate_size_at_first_maturity`) that calculates the "size at the first maturity" because its a function of the mean size of the individuals. For that, we use a linear relationship between mean size and size at the first maturity before/after introducing fishing protection.

      # Stage 2 Trochophore
      stage = 2
      prev_stage = 1

      # immigration
      dispersal_probs1 = (exp(-α[stage, species]) ./ distance_matrix[:, site])
      dispersal_probs1[site] = 0.0
      dispersal_probs1 = dispersal_probs1 .* mig_probs[:, site] .* reduce_trochophore_dispersals
      # dispersal_probs1 = dispersal_probs1 ./ sum(dispersal_probs1)

      # emigration
      dispersal_probs2 = (exp(-α[stage, species]) ./ distance_matrix[site, :])
      dispersal_probs2[site] = 0.0
      dispersal_probs2 = dispersal_probs2 .* mig_probs[site, :] .* reduce_trochophore_dispersals
      # dispersal_probs2 = dispersal_probs2 ./ sum(dispersal_probs2)

      du[site, stage, species] = (r[stage, species] * u[site, prev_stage, species]) -
                        (r[stage+1, species] * u[site, stage, species]) -
                        (d[stage, species] * u[site, stage, species]) +
        (sum(dispersal_probs1 .* u[:, stage, species])) -
        (u[site, stage, species] * sum(dispersal_probs2))

      #stage 3 Veliger
      stage = 3
      prev_stage = 2

      du[site, stage, species] = (r[stage, species] * u[site, prev_stage, species] * ((K[species] - u[site, stage, species]) / K[species])) -
        (r[stage+1, species] * u[site, stage, species]) -
        (d[stage, species] * u[site, stage, species])

      # stage 4 Juvenile
      stage = 4
      prev_stage = 3
      du[site, stage, species] = (r[stage, species] * u[site, prev_stage, species] * ((K[species] - u[site, stage, species]) / K[species])) -
        (r[stage+1, species] * u[site, stage, species]) -
        (d[stage, species] * u[site, stage, species])

      # stage 5 adult
      other_species = species == 1 ? 2 : 1
      stage = 5
      prev_stage = 4
      du[site, stage, species] = (r[stage, species] * u[site, prev_stage, species] * ((K[species] - u[site, stage, species]) / K[species])) -
        (((1 - reproductive_cycle(t)) * exploitation_rates[site, species]) * u[site, stage, species]) -
        (d[stage, species] * u[site, stage, species]) - 
        (competition_coefficient[other_species] * u[site, stage, species] * u[site, stage, other_species])

      # adult sizes
      stage = 6
      du[site, stage, species] = size_growth_rate[species] * u[site, nlifestages+1, species] * (1 - (u[site, nlifestages+1, species] / (size_max[species] - (size_max[species] * ((1 - reproductive_cycle(t)) * exploitation_rates[site, species])))))
    end
  end
  # for i in 1:nsites
  #   for j in 1:(nlifestages-1)
  #     du[i, j] = max(du[i, j], 0.0)
  #   end
  # end
end
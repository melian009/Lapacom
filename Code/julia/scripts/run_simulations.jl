include("Lapacom.jl")
using Main.Lapacom
using Agents

# Define model parameters
nsites = 5
numagents = 100
nspecies = 2
r = 0.4  # reproduction rate
death_probability = 0.04
distMat = rand(nsites,nsites)  # ordered as from row to column
exposure = range(0.9, length=nsites, stop=0.1)
K = rand(20:100, nsites)
tprob = [0.5 for i in 1:5]

# data collection functions
function count_site_pop(model::AbstractModel)
  out = [0 for i in 1:model.nsites]
  for ag in model.agents
    out[ag.pos] += 1
  end
  return Tuple(out)
end

function nspecies1(model)
  counter = 0
  for ag in model.agents
    if ag.species == 1
      counter += 1
    end
  end
  return counter
end

function nspecies2(model)
  counter = 0
  for ag in model.agents
    if ag.species == 2
      counter += 1
    end
  end
  return counter
end

# model running parameters
nsteps = 100
steps_to_collect_data = collect(1:10:nsteps);
propagg = Dict(:model=>[nagents, nspecies1, nspecies2, count_site_pop])

# start the model
model = instantiate_model(distMat=distMat, exposure=exposure, K=K, tprob=tprob, numagents=numagents, nsites=nsites, nspecies=nspecies, death_probability=death_probability, r=r);

# run the model
data = step!(dummystep, model_step!, model, nsteps, propagg, steps_to_collect_data)
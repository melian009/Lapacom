"Function to instantiate the model."
function instantiate_model(;distMat, exposure, K, tprob, death_probability, numagents=30, nsites=7, nspecies=2)
  # An array of arrays for each node of the space.
  agent_positions = [Int32[] for i in 1:nsites]
  # Instantiate the grid structure.
  mygrid = World((1, nsites), grid((1, nsites)), agent_positions)
  # Create a list of agents, each with position (1,1) and one unit of
  # wealth.
  agents = Array{Agent, 1}(undef, numagents)
  for ag in 1:numagents
    agents[ag] = Agent(UInt32(ag), UInt16(1), rand(UInt8(1):UInt8(nspecies)), rand(UInt8(1):UInt8(5)), rand([true, false]))
  end

  # Instantiate and return the model.
  model = Model(mygrid, typeof(agents[1])[], distMat, exposure, K, tprob, numagents, nsites, death_probability)

  # Use the `add_agent!` function to add agents to the model.
  for agent in agents
    add_agent!(agent, model)
  end

  return model
end

"""
Creates a new agent from the two.
"""
function reproduce!(agent1::AbstractAgent, agent2::AbstractAgent, model::AbstractModel)
  if agent1.ismale == agent2.ismale
    throw("Invalid agents selected for mating. They are the same sex.")
  elseif agent1.species != agent2.species
    throw("Invalid agents selected for mating. They are not the same species.")
  end

  newAgent = Agent(typeof(agent1.id)(model.lastId+1), agent1.pos, agent1.species, typeof(agent1.stage)(1), rand([true, false]))
  add_agent!(newAgent, agent1.pos, model)
end

function disperse!(agent::AbstractAgent, model::AbstractModel)
  if agent.stage < 4  # no dispersal in the last two life stages
    # TODO: stage 1 and 2 can disperse globally, stage 3 can disperse no more than regionally.
    new_loc = sample(1:model.nsites, Weights(model.distMat[agent.pos, :]))  # choose a random site given dispersal probabilities.
    move_agent!(agent, new_loc, model)
  end
end

"""
Kills an individual with a probability that is independent of age/site and dependent on the exposure of its habitat to humans. (TODO: correct?)
"""
function death!(model::AbstractModel)
  numagents = nagents(model)
  # TODO: connect death to K.
  death_probs = Weights([model.exposure[i.pos] for i in model.agents])
  todienum = round(Int64, model.deathp * numagents)
  dies_inds = sample(1:numagents, death_probs, todienum, replace=false, ordered=true)
  for index in dies_inds[end:-1:1]
    kill_agent!(model.agents[index], model)
  end
  # TODO: sometimes, it cannot find the agent in model.space.agent_positions[agent.pos]
end

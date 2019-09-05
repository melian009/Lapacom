"Function to instantiate the model."
function instantiate_model(;distMat, exposure, K, tprob, numagents=30, nsites=7, nspecies=2)
  # An array of arrays for each node of the space.
  agent_positions = [Int32[] for i in 1:nsites]
  # Instantiate the grid structure.
  mygrid = World((1, nsites), grid((1, nsites)), agent_positions)
  # Create a list of agents, each with position (1,1) and one unit of
  # wealth.
  agents = Array{Agent, 1}(undef, numagents)
  for ag in 1:numagents
    agents[ag] = Agent(UInt32(ag), (UInt16(1), UInt16(1)), rand(UInt8(1):UInt8(nspecies)), rand(UInt8(1):UInt8(5)), rand([true, false]))
  end

  # Instantiate and return the model.
  model = Model(mygrid, Agent[], distMat, exposure, K, tprob, numagents)

  # Use the `add_agent!` function to add agents to the model.
  for agent in agents
    add_agent!(agent, model)
  end

  return model
end

"""
Creates a new agent from the two.
"""
function reproduce(agent1::Agent, agent2::Agent, model::Model)
  if agent1.ismale == agent2.ismale
    throw("Invalid agents selected for mating. They are the same sex.")
  end

  newAgent = Agent(model.lastId+1, agent1.pos, agent1.species, typeof(agent1.stage)(1), rand([true, false]))
end

function disperse(agent::Agent, model::Model)
end

"""
Kills an individual with a probability that depends on its age and the exposure of its living place to humans.(TODO: correct?)
"""
function death(agent::Agent, model::Model)
end

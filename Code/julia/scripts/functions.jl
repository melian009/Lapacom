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
  model = Model(mygrid, agents, distMat, exposure, K, tprob)

  # Use the `add_agent!` function to add agents to the model.
  for agent in agents
    add_agent!(agent, model)
  end

  return model
end


"Function to instantiate the model."
function instantiate_model(;distMat, exposure, K, tprob, death_probability, r, numagents=30, nsites=7, nspecies=2)
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
  model = Model(mygrid, typeof(agents[1])[], distMat, exposure, K, tprob, numagents, nsites, death_probability, r, nspecies)

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

"""
Pairs males and females within a site. If one of them has a chance more than `model.r`, they will reproduce.
"""
function reproduce!(model::AbstractModel)
  agentnum = nagents(model)
  allprobs = rand(agentnum);
  males, females = separate_sexes(model);
  for sp in 1:model.nspecies
    for siteindex in model.nsites
      if length(males[siteindex][sp]) == 0 || length(females[siteindex][sp]) == 0
        continue
      end
      matches = randmate(males[siteindex][sp], females[siteindex][sp])
      for pair in matches
        if allprobs[pair[1]] <= model.r || allprobs[pair[2]] <= model.r
          agent1 = model.agents[pair[1]]
          agent2 = model.agents[pair[2]]
          reproduce!(agent1, agent2, model)
        end
      end
    end
  end
end

"""
Create two lists, one for males and one for females. Each list has `model.nsites` inner lists which hold agent indices.
"""
function separate_sexes(model::AbstractModel)
  males = [[eltype(model.agents[1].id)[] for j in 1:model.nspecies] for i in 1:model.nsites]
  females = [[eltype(model.agents[1].id)[] for j in 1:model.nspecies] for i in 1:model.nsites]
  for sp in 1:model.nspecies
    for (index, agent) in enumerate(model.agents)
      if agent.ismale && agent.species == sp
        push!(males[agent.pos][sp], index)
      elseif !agent.ismale && agent.species == sp
        push!(females[agent.pos][sp], index)
      end
    end
  end
  return males, females
end

"""
Returns a zip object of n randomly matched pairs between two lists, where n is the size of the smaller list.
"""
function randmate(list1::AbstractArray, list2::AbstractArray)
  jj = eltype(list1)
  list1len = length(list1)
  list2len = length(list2)
  if list1len < list2len
    smallerlist = list1
    largerlist = list2
    smallerlen = list1len
    largerlen = list2len
  else
    smallerlist = list2
    largerlist = list1
    smallerlen = list2len
    largerlen = list1len
  end

  matches = zip(smallerlist[Random.randperm(jj(smallerlen))], largerlist[Random.randperm(jj(largerlen))[1:smallerlen]])
  
  return matches
end

function disperse!(agent::AbstractAgent, model::AbstractModel)
  if agent.stage < 4  # no dispersal in the last two life stages
    new_loc = sample(1:model.nsites, Weights(model.distMat[agent.pos, :]))  # choose a random site given dispersal probabilities.
    move_agent!(agent, new_loc, model)
  end
end

function disperse!(model::AbstractModel)
  for agent in model.agents
    disperse!(agent, model)
  end
end

"""
Kills an individual with a probability that is independent of age/site and dependent on the exposure of its habitat to humans.
"""
function kill!(model::AbstractModel)
  numagents = nagents(model)
  # TODO: connect death to K for later versions
  death_probs = Weights([model.exposure[i.pos] for i in model.agents])
  todienum = round(Int64, model.deathp * numagents)
  dies_inds = sample(1:numagents, death_probs, todienum, replace=false, ordered=true)
  for index in dies_inds[end:-1:1]
    kill_agent!(model.agents[index], model)
  end
  # FIXME: there are mismatches between agent.pos and model.space.agent_positions
end


function model_step!(model::AbstractModel)
  reproduce!(model)
  disperse!(model)
  kill!(model)
end
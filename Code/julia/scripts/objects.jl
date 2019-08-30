
mutable struct Agent{T<:Integer, U<:Integer, V<:Integer} <: AbstractAgent
  id::T
  pos::V
  species::U  # species ID
  stage::U  # life stage ID
end

mutable struct World{T<:Integer, V<:AbstractArray, W<:SimpleGraph} <: AbstractSpace
  agent_positions::V
  space::W
  dimensions::T
end

mutable struct Model{T<:AbstractFloat, U<:AbstractSpace, V<:AbstractArray, W<:AbstractArray} <: AbstractModel
  space::U
  agents::V
  distMat::Array{T, 2}  # distance between each pair of patches
  exposure::Array{T, 1}  # an array of amount of exposure for each site
  K::W  # carrying capacity at each site
end

"Function to instantiate the model."
function instantiate_model(; numagents=30, griddims=(1,7))
  # An array of arrays for each node of the space.
  agent_positions = [Int32[] for i in 1:gridsize(griddims)]
  # Instantiate the grid structure.
  mygrid = World(griddims, grid(griddims), agent_positions)
  # Create a list of agents, each with position (1,1) and one unit of
  # wealth.
  agents = [Agent(i, (1,1), 1) for i in 1:numagents]

  # Instantiate and return the model.
  model = Model(mygrid, MyAgent[], random_activation)

  # Use the `add_agent!` function to add agents to the model.
  for agent in agents
    add_agent!(agent, model)
  end

  return model
end
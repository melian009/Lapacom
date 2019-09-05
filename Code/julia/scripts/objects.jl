
"An object holding agent parameters."
mutable struct Agent{T<:Integer, U<:Integer, V<:Integer} <: AbstractAgent
  id::T
  pos::V
  species::U  # species ID
  stage::U  # life stage ID
  ismale::Bool
end

"An object describing the world of simulation"
mutable struct World{T<:Tuple{Integer, Integer}, V<:AbstractArray, W<:SimpleGraph} <: AbstractSpace
  dimensions::T
  space::W
  agent_positions::V
end

"An object holding model parameters"
mutable struct Model{T<:AbstractArray, U<:AbstractSpace, V<:AbstractArray, W<:AbstractArray, X<:AbstractArray, Y<:AbstractArray, Z<:Integer, A<:Integer, B<:AbstractFloat} <: AbstractModel
  space::U
  agents::V
  distMat::T  # distance between each pair of patches. ordered as from row to column
  exposure::Y  # an array of amount of exposure for each site
  K::W  # carrying capacity at each site
  tprob::X  # transition probability between angents' life stages.
  lastId::Z  # last agent ID used in the model.
  nsites::A  # number of sites
  deathp::B  # global death probability of agents, which is irrespective of age/size and location
end
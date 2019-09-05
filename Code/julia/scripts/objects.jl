
"An object holding agent parameters."
mutable struct Agent{T<:Integer, U<:Integer, V<:Tuple{Integer, Integer}} <: AbstractAgent
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
mutable struct Model{T<:AbstractArray, U<:AbstractSpace, V<:AbstractArray, W<:AbstractArray, X<:AbstractArray, Y<:AbstractArray, Z<:Integer} <: AbstractModel
  space::U
  agents::V
  distMat::T  # distance between each pair of patches
  exposure::Y  # an array of amount of exposure for each site
  K::W  # carrying capacity at each site
  tprob::X  # transition probability between angents' life stages.
  lastId::Z  # last agent ID used in the model.
end
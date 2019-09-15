module Lapacom

using Agents
using DrWatson
using StatsBase: sample, Weights
include("objects.jl")
include("functions.jl")

export Agents, Agent, World, Model, instantiate_model, model_step!, kill!, disperse!, reproduce!

end  # module
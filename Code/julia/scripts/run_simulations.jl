include("Lapacom.jl")
using Main.Lapacom
using Agents

# Define model parameters
nsites = 5
numagents = 100
nspecies = 2
r = 0.4  # reproduction rate
death_probability = 0.02
distMat = rand(nsites,nsites)  # ordered as from row to column
exposure = range(0.9, length=nsites, stop=0.1)
K = rand(20:100, nsites)
tprob = [0.5 for i in 1:5]

# start the model
model = instantiate_model(distMat=distMat, exposure=exposure, K=K, tprob=tprob, numagents=numagents, nsites=nsites, nspecies=nspecies, death_probability=death_probability, r=r);

# run the model

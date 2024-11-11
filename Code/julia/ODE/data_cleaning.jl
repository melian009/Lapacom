using Pkg
Pkg.activate(".")
using DataFrames
using CSV
using Dates
include("distancematrix.jl")

uncleaned_file = "../../../Data/ToAnalyze/Madeira/BD_LIMPETS_MAD_1996-2018.csv"

df = CSV.read(uncleaned_file, DataFrame)

desired_columns = ["Individual_number", "Species", "Year", "Month", "Total_length_mm", "Weight_g", "Sampling_site", "Lat", "Long", "Protective_regime", "Proximity_human_settlements"]

df = df[.!ismissing.(df.Lat), desired_columns]
rename!(lowercase, df)
rename!(df, Dict(:long => :lon, :individual_number => :id))

df.date = Date.(df.year, df.month)

# Clean site names
df.sampling_site = strip.(df.sampling_site)
df[df.sampling_site .== "Faj\x8b dos Padres", :sampling_site] .= "Faja Padres"
desired_sites = ["Porto Moniz", "Pa\x9cl do Mar", "S\x8bo Vicente", "Ribeira Brava", "Santa Cruz", "Funchal", "Cani\x8dal", "Desertas"]
df = df[findall(x -> in(x, desired_sites), df.sampling_site), :]

@assert length(unique(df.id)) == size(df, 1) "There are non-unique id"
@assert !(any(x -> ismissing(x), df.sampling_site)) "There are missing sampling sites"

sort!(df, [:date, :total_length_mm])

df.lat = dms2decimal.(split_dms_str.(df.lat))
df.lon = dms2decimal.(split_dms_str.(df.lon))

CSV.write("data.csv", df)

df = CSV.read("data.csv", DataFrame)

## Create distance matrix
distance_matrix = create_distance_matrix("data.csv"; loc_colname = :sampling_site, lat_colname = :lat, lon_colname = :lon)
CSV.write("distance_matrix.csv", distance_matrix)
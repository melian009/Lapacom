using Pkg
Pkg.activate(".")
using DataFrames
using CSV
using Dates

uncleaned_file = "../../../Data/ToAnalyze/Madeira/BD_LIMPETS_MAD_1996-2018.csv"

df = CSV.read(uncleaned_file, DataFrame)

desired_columns = ["Individual_number", "Species", "Year", "Month", "Total_length_mm", "Weight_g", "Sampling_site", "Lat", "Long", "Protective_regime", "Proximity_human_settlements"]

df = df[.!ismissing.(df.Lat), desired_columns]
rename!(lowercase, df)
rename!(df, Dict(:long => :lon, :individual_number => :id))

df.date = Date.(df.year, df.month)

@assert length(unique(df.id)) == size(df, 1) "There are non-unique id"
@assert !(any(x -> ismissing(x), df.sampling_site)) "There are missing sampling sites"

sort!(df, [:date, :total_length_mm])

"""
    dms2decimal(degrees::Int, minutes::AbstractFloat, seconds::AbstractFloat, direction::AbstractString)

Convert geographic coordinates from DMS (degrees, minutes, seconds) to decimal coordinates.

## Formula:

* DEC = (DEG + (MIN * 1/60) + (SEC * 1/60 * 1/60))
* Assumes S/W are negative. 

## Source

https://gist.github.com/jeteon/89c41e4081d87b798d8006b16a52c695
"""
function dms2decimal(degrees::Int, minutes::AbstractFloat, seconds::AbstractFloat, direction::AbstractString)
  dsign = direction in ["S", "s", "W", "w"] ? -1 : 1
  DEC = dsign * (degrees + (minutes * 1 / 60) + (seconds * 1 / 60 * 1 / 60))
end

dms2decimal(x::Tuple) = dms2decimal(x[1], x[2], x[3], x[4])

"""
This function is written specifically for the DMS format in the `df`.
"""
function split_dms_str(dms_str::AbstractString)
  o = eachmatch(r"[\d|\.|(W|w|S|s|N|n|E|e)]+", dms_str)
  oo = collect(o)
  return parse(Int, oo[1].match), parse(Float64, oo[2].match), parse(Float64, oo[3].match), oo[4].match
end

df.lat = dms2decimal.(split_dms_str.(df.lat))
df.lon = dms2decimal.(split_dms_str.(df.lon))

CSV.write("data.csv", df)

df = CSV.read("data.csv", DataFrame)

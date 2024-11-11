function haversine(lat1, lon1, lat2, lon2)
  π180 = π / 180.0
  # distance between latitudes and longitudes 
  dLat = (lat2 - lat1) * π180
  dLon = (lon2 - lon1) * π180

  # convert to radians 
  lat1 *= π180
  lat2 *= π180

  # apply formulae 
  a = (sin(dLat / 2)^2) +
      (sin(dLon / 2)^2) *
      cos(lat1) * cos(lat2)
  rad = 6371
  c = 2 * asin(sqrt(a))
  return rad * c
end

function create_distance_matrix(inputfile; loc_colname::Symbol=:sampling_site, lat_colname::Symbol=:lat, lon_colname::Symbol=:lon)
  df = CSV.read(inputfile, DataFrame)
  df = df[:, [loc_colname, lat_colname, lon_colname]]
  rows = [findfirst(x -> x == i, df[:, loc_colname]) for i in unique(df[:, loc_colname])]
  df = df[rows, :]

  dist_mat = zeros(size(df, 1), size(df, 1))
  for r1 in 1:size(df, 1)
    for r2 in 1:size(df, 1)
      dist = haversine(df[r1, lat_colname], df[r1, lon_colname], df[r2, lat_colname], df[r2, lon_colname])
      dist_mat[r1, r2] = dist
    end
  end
  colnames = df[:, loc_colname]
  dist_df = DataFrame(dist_mat, :auto)
  rename!(dist_df, Symbol.(colnames))
  dist_df[:, :site] = colnames
  return dist_df
end

"""
    dms2decimal(degrees::Int, minutes::AbstractFloat, seconds::AbstractFloat, direction::AbstractString)

Convert geographic coordinates from DMS (degrees, minutes, seconds) to decimal coordinates.

## Formula:

* decimal = degrees + (minutes * 1/60) + (seconds * 1/3600)
* Assumes S/W are negative. 

## Source

https://www.rapidtables.com/convert/number/degrees-minutes-seconds-to-degrees.html
"""
function dms2decimal(degrees::Int, minutes::AbstractFloat, seconds::AbstractFloat, direction::AbstractString)
  dsign = direction in ["S", "s", "W", "w"] ? -1 : 1
  DEC = dsign * (degrees + (minutes * 1 / 60) + (seconds * 1 / 3600))
end

dms2decimal(x::Tuple) = dms2decimal(x[1], x[2], x[3], x[4])


"""
Split a string containing the coordinates in DMS to separate numbers.
"""
function split_dms_str(dms_str::AbstractString)
  o = eachmatch(r"[\d|\.|(W|w|S|s|N|n|E|e)]+", dms_str)
  oo = collect(o)
  return parse(Int, oo[1].match), parse(Float64, oo[2].match), parse(Float64, oo[3].match), oo[4].match
end

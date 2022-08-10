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

function create_distance_matrix(inputfile; loc_colname::Symbol = :sampling_site, lat_colname::Symbol= :lat, lon_colname::Symbol = :lon)
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


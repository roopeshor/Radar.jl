"Read a tab-separated radiosonde table and retain the ascending, unique-altitude segment."
function read_radiosonde(filepath::AbstractString)
  isfile(filepath) ||
    throw(ArgumentError("Radiosonde file not found at: $filepath"))
  lines = readlines(filepath)
  isempty(lines) && throw(ArgumentError("Radiosonde file is empty: $filepath"))
  rows = Vector{Vector{Float64}}()
  for line in lines[2:end] # first row is the table header
    isempty(strip(line)) && continue
    fields = split(chomp(line), '\t')
    length(fields) >= 10 || continue
    vals = [something(tryparse(Float64, strip(fields[i])), NaN) for i = 1:10]
    push!(rows, vals)
  end
  isempty(rows) && throw(
    ArgumentError("No numeric radiosonde observations found in $filepath"),
  )
  raw = reduce(vcat, permutedims.(rows))
  raw_time = raw[:, 1]
  raw_p = raw[:, 3]
  raw_t = raw[:, 4]
  raw_ws = raw[:, 6]
  raw_wd = raw[:, 7]
  raw_u = raw[:, 8]
  raw_v = raw[:, 9]
  raw_alt = raw[:, 10]
  valid =
    isfinite.(raw_alt) .& isfinite.(raw_u) .& isfinite.(raw_v) .& (raw_alt .> 0)
  valid_indices = findall(valid)
  isempty(valid_indices) &&
    throw(ArgumentError("No valid positive-altitude radiosonde points"))
  max_alt_idx = argmax(ifelse.(isfinite.(raw_alt), raw_alt, -Inf))
  ascent = [i for i in valid_indices if i <= max_alt_idx]

  # Preserve the first occurrence of each altitude, as MATLAB unique(...,'stable') does.
  seen = Set{Float64}()
  unique_indices = Int[]
  for i in ascent
    if !(raw_alt[i] in seen)
      push!(seen, raw_alt[i])
      push!(unique_indices, i)
    end
  end
  order = sortperm(raw_alt[unique_indices])
  idx = unique_indices[order]
  altitude = raw_alt[idx]
  @info "Loaded radiosonde data" points=length(altitude) min_altitude=first(
    altitude,
  ) max_altitude=last(altitude)
  return (;
    altitude_m = altitude,
    u = raw_u[idx],
    v = raw_v[idx],
    ws = raw_ws[idx],
    wd = raw_wd[idx],
    temp = raw_t[idx],
    time = raw_time[idx],
    pressure = raw_p[idx],
    filepath = String(filepath),
  )
end

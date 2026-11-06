function _read_numeric_table(path::AbstractString)
  lines = readlines(path)
  isempty(lines) && throw(ArgumentError("Empty table: $path"))
  delim = occursin('\t', first(lines)) ? '\t' : ','
  headers = strip.(split(chomp(first(lines)), delim))
  rows = Vector{Vector{Float64}}()
  for line in lines[2:end]
    isempty(strip(line)) && continue
    fields = strip.(split(chomp(line), delim))
    length(fields) < length(headers) && continue
    vals = Float64[]
    for f in fields[1:length(headers)]
      value = tryparse(Float64, replace(f, '"' => ""))
      push!(vals, value === nothing ? NaN : value)
    end
    push!(rows, vals)
  end
  data =
    isempty(rows) ? zeros(Float64, 0, length(headers)) :
    reduce(vcat, permutedims.(rows))
  return headers, data
end

function _column(headers, data, name::AbstractString)
  i = findfirst(==(name), headers)
  i === nothing && throw(
    ArgumentError(
      "Required column '$name' not found. Available columns: $(join(headers, ", "))",
    ),
  )
  return vec(data[:, i])
end

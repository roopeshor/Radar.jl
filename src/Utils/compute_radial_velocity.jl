function compute_radial_velocity(
  u::AbstractVector,
  v::AbstractVector,
  w::AbstractVector,
  az::Real,
  oz::Real,
)
  length(u) == length(v) == length(w) ||
    throw(DimensionMismatch("u, v and w must have equal lengths"))
  if iszero(oz)
    return Float64.(w)
  end
  azr = deg2rad(az)
  ozr = deg2rad(oz)
  return (u .* sin(azr) .+ v .* cos(azr)) .* sin(ozr) .+ w .* cos(ozr)
end

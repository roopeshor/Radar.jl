function compute_height_bins(
  start_range::Real,
  end_range::Real,
  num_range_bins::Integer,
)
  num_range_bins >= 2 ||
    throw(ArgumentError("num_range_bins must be at least 2"))
  return collect(
    range(Float64(start_range), Float64(end_range), length = num_range_bins),
  )
end

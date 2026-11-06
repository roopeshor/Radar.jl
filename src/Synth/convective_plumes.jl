function convective_plumes(
  h::AbstractVector;
  plume_count::Integer = 4,
  plume_max_alt::Real = 6.0,
)
  z = Float64.(h)
  out = zeros(length(z))
  for _ = 1:plume_count
    center = 0.8 + rand()*(plume_max_alt-0.8)
    peak = (rand()-0.3)*2.5
    width = 0.2 + rand()*0.4
    out .+= peak .* exp.(-((z .- center) ./ width) .^ 2)
  end
  return out
end

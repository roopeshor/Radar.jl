function somali_jet(
  h::AbstractVector;
  sj_w::Real = 12.0,
  sj_uv::Real = 14.0,
  sj_center::Real = 1.5,
  sj_s::Real = 0.8,
)
  shape = exp.(-((h .- sj_center) ./ sj_s) .^ 2)
  return (sj_w .* shape, sj_uv .* shape)
end

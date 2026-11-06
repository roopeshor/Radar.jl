function tropical_easterly_jet(
  h::AbstractVector;
  tej_u::Real = 10.0,
  tej_v::Real = 3.0,
  tej_alt::Real = 15.0,
  tej_vthick::Real = 2.5,
)
  shape = exp.(-((h .- tej_alt) ./ tej_vthick) .^ 2)
  return (tej_u .* shape, tej_v .* shape)
end

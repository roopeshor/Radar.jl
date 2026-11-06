function orographic_updraft(
  h::AbstractVector,
  u_mean::AbstractVector;
  ou_alpha_slope::Real = atan(3/20),
  ou_H_lift::Real = 2.0,
)
  return max.(0.0, abs.(u_mean)) .* tan(ou_alpha_slope) .* exp.(-h ./ ou_H_lift)
end

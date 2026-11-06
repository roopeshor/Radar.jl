function synoptic_bg(
  h::AbstractVector;
  bg_z0_m::Real = 0.05,
  bg_ustar::Real = 0.3,
  bg_kappa::Real = 0.4,
  bg_blh::Real = 1.5,
)
  z_m = h .* 1000
  return (bg_ustar/bg_kappa) .* log.((z_m .+ bg_z0_m) ./ bg_z0_m) .*
         exp.(-h ./ bg_blh)
end

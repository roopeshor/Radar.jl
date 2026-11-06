function gravity_waves(
  h::AbstractVector;
  gw_nmodes::Integer = 20,
  gw_amp::Real = 2.5,
  gw_max_lambda::Real = 4.0,
  gw_damp_scale::Real = 12.0,
  gw_cf::Real = 0.05,
  lw_Hs::Real = 7.0,
)
  z = Float64.(h)
  n = length(z)
  m = collect(range(2π/gw_max_lambda, 2π/0.4, length = gw_nmodes))
  phase_u = 2π .* rand(gw_nmodes)
  phase_v = 2π .* rand(gw_nmodes)
  energy = m .^ (-1.5)
  energy ./= sum(energy)
  uw = zeros(n)
  vw = zeros(n)
  for i = 1:gw_nmodes
    uw .+= sqrt(energy[i]) .* cos.(m[i] .* z .+ phase_u[i])
    vw .+= sqrt(energy[i]) .* sin.(m[i] .* z .+ phase_v[i])
  end
  amp = exp.(z ./ (2lw_Hs)) .* exp.(-z ./ gw_damp_scale)
  su = std(uw)
  sv = std(vw)
  su > 0 && (uw ./= su)
  sv > 0 && (vw ./= sv)
  u = gw_amp .* uw .* amp
  v = gw_amp .* vw .* amp
  return (u, v, -gw_cf .* u)
end

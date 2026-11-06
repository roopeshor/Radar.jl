function compute_spectrum_broadening(
  z::AbstractVector,
  u::AbstractVector,
  v::AbstractVector,
  w::AbstractVector,
  v_radial::AbstractVector;
  sigma_turb_base::Real = 0.2,
  sigma_chi::Real = 0.2,
)
  n = length(z)
  n >= 2 || throw(ArgumentError("At least two altitude gates are required"))
  dz_m = (z[2]-z[1])*1000
  v_cross = sqrt.(max.(u .^ 2 .+ v .^ 2 .+ w .^ 2 .- v_radial .^ 2, 0))
  dv_dz = abs.(gradient(v_radial, dz_m))
  sigma_shear = dv_dz .* abs(dz_m) ./ sqrt(12)
  theta = deg2rad(3.3)
  sigma_beam = theta/(2sqrt(2log(2))) .* abs.(v_cross)
  sigma_base =
    sqrt.(sigma_shear .^ 2 .+ sigma_beam .^ 2 .+ Float64(sigma_turb_base)^2)
  μ = -0.5sigma_chi^2
  noise = randn(n)
  corr_gates = max(3, round(Int, 0.45/(z[2]-z[1])))
  eta = gaussian_smooth(noise, corr_gates)
  s = std(eta)
  s > 0 && (eta ./= s)
  return sigma_base .* exp.(μ .+ sigma_chi .* eta)
end

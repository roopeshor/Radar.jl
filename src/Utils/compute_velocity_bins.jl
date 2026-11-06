function compute_velocity_bins(nfft::Integer, ipp::Real, n_coh::Real)
  nfft > 0 || throw(ArgumentError("nfft must be positive"))
  c = 299792458.0
  fc = 205e6
  λ = c / fc
  effective_PRI = ipp * 1e-6 * n_coh
  fd = collect((-fld(nfft, 2)):(cld(nfft, 2) - 1)) ./ (nfft * effective_PRI)
  return λ .* fd ./ 2
end

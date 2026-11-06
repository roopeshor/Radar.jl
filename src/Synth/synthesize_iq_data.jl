function synthesize_iq_data(;
  vr::AbstractVector,
  Header::RadarData,
  spec_w::AbstractVector,
  SNR::AbstractVector = Float64[],
  add_iq_noise::Bool = true,
  dc_rel_power::Real = 0.0,
  gc_rel_power::Real = 0.005,
  sigma_dc::AbstractVector = Float64[],
  sigma_gc::AbstractVector = Float64[],
  bright_bands::AbstractMatrix = zeros(0, 2),
)
  nfft = Int(Header.nFFT)
  nr = Int(Header.nRangeBins)
  ipp = Header.ipp_us
  coh = Header.n_coh
  length(vr) == nr ||
    throw(DimensionMismatch("vr length must match header range bins"))
  length(spec_w) == nr ||
    throw(DimensionMismatch("spec_w length must match header range bins"))
  pnoise = 1000.0
  wavelength = 3e8/205e6
  dbs_v = wavelength/2
  fd = -Float64.(vr) ./ dbs_v
  prf = 1/(ipp*1e-6*coh)
  sigma_f = 2 .* Float64.(spec_w) ./ wavelength
  any(sigma_f .<= 0) && (sigma_f .= max.(sigma_f, eps(Float64)))
  isempty(SNR) && (SNR = 10 .^ (collect(range(40.0, -40.0, length = nr)) ./ 10))
  isempty(sigma_dc) && (sigma_dc = ones(nr) ./ 20)
  isempty(sigma_gc) && (sigma_gc = 1 ./ (collect(1:nr) ./ 10) ./ 5)
  f = collect(range(-prf/2, prf/2-prf/nfft, length = nfft))
  frow = reshape(f, 1, nfft)
  fdcol = reshape(fd, nr, 1)
  sigcol = reshape(sigma_f, nr, 1)
  snrcol = reshape(Float64.(SNR), nr, 1)
  sdc = reshape(Float64.(sigma_dc), nr, 1)
  sgc = reshape(Float64.(sigma_gc), nr, 1)
  psig = pnoise .* snrcol
  ssig = exp.(-((frow .- fdcol) .^ 2) ./ (2 .* sigcol .^ 2))
  sums = sum(ssig; dims = 2)
  scale = zeros(nr, 1)
  valid = sums .> 0
  scale[valid] .= (psig[valid] ./ sums[valid]) ./ (sqrt(2π) .* sigcol[valid])
  ssig .*= scale
  spectrum = ssig .+ pnoise/nfft
  row_sums = sum(spectrum; dims = 2)
  spectrum .+= (dc_rel_power .* row_sums) .* exp.(-frow .^ 2 ./ (2 .* sdc .^ 2))
  row_sums = sum(spectrum; dims = 2)
  spectrum .+= (gc_rel_power .* row_sums) .* exp.(-frow .^ 2 ./ (2 .* sgc .^ 2))

  for k in axes(bright_bands, 1)
    idx = round(Int, bright_bands[k, 1])
    if 1 <= idx <= nr
      amp_scale = -3.5 + (2.0+3.5)*(k-1)/max(nr-1, 1)
      amp = sum(@view spectrum[idx, :]) * 10^amp_scale
      spectrum[idx, :] .+= amp .* exp.(-f .^ 2 ./ (2*bright_bands[k, 2]^2))
    end
  end
  sroot = sqrt.(spectrum ./ 2)
  z = if add_iq_noise
    sroot .* (randn(nr, nfft) .+ im .* randn(nr, nfft))
  else
    ComplexF64.(sroot)
  end
  time_data = ifft(ifftshift(z, 2), 2) .* sqrt(nfft)
  rounded =
    complex.(_round_away.(real.(time_data)), _round_away.(imag.(time_data)))
  return reshape(ComplexF64.(rounded), nr, nfft, 1)
end

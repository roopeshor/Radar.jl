function HS_noise_estimate(beam::RadarData)
  # Faithful translation of the MATLAB implementation. Note that its P_n
  # denominator starts at N (rather than 1); this is unusual, but retained
  # so MATLAB/Julia comparisons use the same estimator.
  power = sort(abs2.(beam.spectra); dims = 2)
  nr, nfft = size(power)
  M = beam.m_sNumOfInCohIntegrations

  p_den = reshape(nfft:(2nfft - 1), 1, nfft)
  p_n = cumsum(power; dims = 2) ./ p_den
  q_den = reshape(1:nfft, 1, nfft)
  q_n = cumsum(power .^ 2; dims = 2) ./ q_den .- p_n .^ 2
  r_n = p_n .^ 2 ./ (q_n .* M)

  noise = zeros(Float64, nr)
  for i = 1:nr
    idx_from_right =
      findfirst(x -> !isnan(x) && x > 1.0, reverse(vec(r_n[i, :])))
    idx0 = idx_from_right === nothing ? 0 : idx_from_right - 1
    valid_idx0 = nfft - 1 - idx0
    noise[i] = p_n[i, valid_idx0 + 1]
  end
  return noise
end

"Calculate spectral moments for one peak after subtracting a scalar noise floor."
function calculate_moments(
  P::AbstractVector,
  vel_axis::AbstractVector,
  noise_level::Real,
)
  length(P) == length(vel_axis) ||
    throw(DimensionMismatch("P and vel_axis lengths differ"))
  ps = max.(Float64.(P) .- noise_level, 0.0)
  maxpow, peak_idx = findmax(ps)
  if maxpow <= 0
    return (M0 = 0.0, M1 = NaN, M2 = NaN, width = NaN, SNR = -Inf)
  end
  left = peak_idx
  while left > 1 && ps[left - 1] > 0
    left -= 1
  end
  right = peak_idx
  while right < length(ps) && ps[right + 1] > 0
    right += 1
  end
  pseg = @view ps[left:right]
  vseg = @view vel_axis[left:right]
  m0 = sum(pseg)
  m0 <= 0 && return (M0 = 0.0, M1 = NaN, M2 = NaN, width = NaN, SNR = -Inf)
  m1 = sum(vseg .* pseg) / m0
  m2 = sum((vseg .- m1) .^ 2 .* pseg) / m0
  total_noise = noise_level * length(ps)
  snr = total_noise > 0 ? 10log10(m0 / total_noise) : Inf
  return (M0 = m0, M1 = m1, M2 = m2, width = sqrt(m2), SNR = snr)
end

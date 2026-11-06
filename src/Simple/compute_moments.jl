function compute_moments(beam::RadarData)
  p =
    isempty(beam.denoised_spectra) ? abs2.(beam.spectra) : beam.denoised_spectra
  nr, nfft = size(p)
  m0 = zeros(Float64, nr)
  m1 = fill(NaN, nr)
  m2 = fill(NaN, nr)
  snr = fill(NaN, nr)
  for r = 1:nr
    row = @view p[r, :]
    peak = argmax(row)
    left = peak
    right = peak
    while left > 1 && row[left - 1] > 0
      left -= 1
    end
    while right < nfft && row[right + 1] > 0
      right += 1
    end
    block = row[left:right]
    m0[r] = sum(block)
    m0[r] <= 0 && continue
    idx = (left:right) .- 1
    freq = (idx .- nfft/2) ./ (beam.ipp_us * beam.n_coh * nfft * 1e-6)
    m1[r] = sum(block .* freq) / m0[r]
    m2[r] = sum((freq .- m1[r]) .^ 2 .* block) / m0[r]
    # Kept faithful to the MATLAB implementation, where SNR stores the third central moment.
    snr[r] = sum((freq .- m1[r]) .^ 3 .* block) / m0[r]
  end
  beam.M0 = m0
  beam.M1 = m1
  beam.M2 = m2
  beam.SNR = snr
  return beam
end

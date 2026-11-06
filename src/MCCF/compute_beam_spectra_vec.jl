function compute_beam_spectra_vec(beam::RadarData; max_axis = 2)
  data = beam.BeamData
  nr, nfft, ni = size(data)
  m = nr * ni
  seg_len = 512
  num_delays = 4
  nfft >= seg_len + num_delays - 1 ||
    throw(ArgumentError("NFFT must be large enough for MCCF delay segments"))
  step = fld(nfft-seg_len, num_delays-1)
  signals = reshape(permutedims(data, (1, 3, 2)), m, nfft)
  f = fft(signals, 2)
  noise = [median(abs.(@view f[r, :])) for r = 1:m]
  mask = [1, 2, 3, nfft, nfft-1]
  for r = 1:m
    f[r, mask] .= noise[r] .* cis.(2π .* rand(length(mask)))
  end
  clean = ifft(f, 2)
  lconv = 2seg_len-1
  first_f = fft_n(clean[:, 1:seg_len], lconv, 2)
  sum_f = zeros(ComplexF64, m, lconv)
  for d = 2:num_delays
    start = 1 + (d-1)*step
    seg = conj.(reverse(clean[:, start:(start + seg_len - 1)]; dims = 2))
    sum_f .+= fft_n(seg, lconv, 2)
  end
  conv_accum = ifft(first_f .* sum_f, 2) ./ (num_delays-1)
  padded = zeros(ComplexF64, m, nfft)
  if lconv < nfft
    off = fld(nfft-lconv, 2)
    padded[:, (off + 1):(off + lconv)] .= conv_accum
  elseif lconv > nfft
    off = fld(lconv-nfft, 2)
    padded .= conv_accum[:, (off + 1):(off + nfft)]
  else
    padded .= conv_accum
  end
  smoothed = gaussian_smooth(padded, 5; dims = 2)
  spec = abs2.(fftshift(fft(smoothed, 2), 2))
  spec3 = reshape(spec, nr, ni, nfft)
  spectra = dropdims(mean(spec3; dims = 2); dims = 2)
  beam.denoised_spectra = _normalise_spectra(spectra, max_axis)
  return beam
end

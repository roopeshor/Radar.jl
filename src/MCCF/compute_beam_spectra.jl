function compute_beam_spectra(beam::RadarData; max_axis = 2)
  data = beam.BeamData
  nr, nfft, ni = size(data)
  seg_len = 512
  num_delays = 4
  nfft >= seg_len + num_delays - 1 ||
    throw(ArgumentError("NFFT must be large enough for MCCF delay segments"))
  step = fld(nfft-seg_len, num_delays-1)
  cube = zeros(Float64, nr, nfft, ni)
  mask = [1, 2, 3, nfft, nfft-1]
  for r = 1:nr, q = 1:ni
    signal = vec(data[r, :, q])
    fsig = fft(signal)
    noise = median(abs.(fsig))
    fsig[mask] .= noise .* cis.(2π .* rand(length(mask)))
    clean = ifft(fsig)
    first_seg = clean[1:seg_len]
    accum = zeros(ComplexF64, 2seg_len-1)
    for d = 2:num_delays
      start = 1 + (d-1)*step
      seg = clean[start:(start + seg_len - 1)]
      accum .+= _conv1d(first_seg, conj.(reverse(seg)))
    end
    accum ./= num_delays-1
    padded = zeros(ComplexF64, nfft)
    lconv = length(accum)
    if lconv < nfft
      off = fld(nfft-lconv, 2)
      padded[(off + 1):(off + lconv)] .= accum
    elseif lconv > nfft
      off = fld(lconv-nfft, 2)
      padded .= accum[(off + 1):(off + nfft)]
    else
      padded .= accum
    end
    smoothed = gaussian_smooth(padded, 5)
    cube[r, :, q] .= abs2.(fftshift(fft(smoothed)))
  end
  spectra = dropdims(mean(cube; dims = 3); dims = 3)
  beam.denoised_spectra = _normalise_spectra(spectra, max_axis)
  return beam
end

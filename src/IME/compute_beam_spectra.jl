function compute_beam_spectra(beam::RadarData)
  raw = beam.BeamData
  isempty(raw) && throw(ArgumentError("BeamData is empty"))
  _, nfft, ni = size(raw)
  window = reshape(hanning_window(nfft), 1, nfft, 1)
  windowed = raw .* window
  transformed = fftshift(fft(windowed, 2), 2)
  beam.denoised_spectra =
    dropdims(mean(abs2.(transformed); dims = 3); dims = 3) ./ ni
  return beam
end

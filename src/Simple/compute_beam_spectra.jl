function compute_beam_spectra(beam::RadarData)
  isempty(beam.BeamData) && throw(ArgumentError("BeamData is empty"))
  centered = beam.BeamData .- mean(beam.BeamData; dims = 2)
  cube = abs2.(fftshift(fft(centered, 2), 2))
  spectra = dropdims(mean(cube; dims = 3); dims = 3)
  peak = maximum(spectra)
  beam.denoised_spectra =
    peak == 0 ? zeros(Float64, size(spectra)) : Float64.(spectra ./ peak)
  return beam
end

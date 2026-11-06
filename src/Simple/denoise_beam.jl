function denoise_beam(beam::RadarData)
  spec =
    isempty(beam.denoised_spectra) ? abs2.(beam.spectra) : beam.denoised_spectra
  noise = isempty(beam.noise_level) ? HS_noise_estimate(beam) : beam.noise_level
  clean = max.(spec .- reshape(noise, :, 1), 0.0)
  beam.noise_level = noise
  beam.denoised_spectra = clean
  return beam
end

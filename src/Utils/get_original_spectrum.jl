function get_original_spectrum(beam::RadarData)
  isempty(beam.BeamData) && return zeros(ComplexF64, 0, 0)
  spectra_cube = fftshift(fft(beam.BeamData, 2), 2)
  spectra = dropdims(mean(spectra_cube; dims = 3); dims = 3)
  peak = maximum(abs, spectra)
  return peak == 0 ? zeros(ComplexF64, size(spectra)) :
         ComplexF64.(spectra ./ peak)
end

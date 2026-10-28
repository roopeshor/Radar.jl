
"""
    get_original_spectrum(beam)

Computes normalized spectrum from given beam with no further processing
doesnt do any squaring

# Arguments
- `beam::RadarData`: Radar beam object containing raw IQ complex time-series cube (BeamData).

# Returns
- spectra (:, 1024) double
"""
function get_original_spectrum(beam::RadarData)
    spectraCube = fft(beam.BeamData, 2)
    spectraCube = fftshift(spectraCube, 2)

    # Incoherent integration and normalization
    # taking mean right now works because there is onty 1 incoherent integration
    spectraCube = mean(spectraCube, dims=3)
    spectra = dropdims(spectraCube, dims=3)

    return spectra ./ maximum(spectra)
end

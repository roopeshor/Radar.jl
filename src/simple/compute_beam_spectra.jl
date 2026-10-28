
"""
    compute_beam_spectra(beam)

Computes Doppler power spectra using mean DC subtraction and FFT.

# Arguments
- `beam::RadarData`: RadarData object

# Returns
- `beam::RadarData`
"""
function compute_beam_spectra(beam::RadarData)
    # DC Removal in time domain
    beamData = beam.BeamData .- mean(beam.BeamData, dims=2)

    # Doppler FFT across pulse time series
    spectraCube = fft(beamData, 2)
    spectraCube = abs.(fftshift(spectraCube, 2)) .^ 2

    # Incoherent integration and global peak normalization
    spectraCube = mean(spectraCube, dims=3)
    spectra = dropdims(spectraCube, dims=3)
    beam.denoised_spectra = spectra ./ maximum(spectra)

    return beam
end

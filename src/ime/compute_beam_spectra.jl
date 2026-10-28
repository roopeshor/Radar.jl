"""
    compute_beam_spectra(beam)

Takes raw complex IQ data and generates the averaged Doppler power spectrum.

# Arguments
- `beam::RadarData`: RadarData beam object.

# Returns
- `beam::RadarData`
"""
function compute_beam_spectra(beam::RadarData)
    raw_iq = beam.BeamData
    _, nfft, num_incoh = size(raw_iq)

    # Hann window (manual — no DSP package required)
    win = reshape(0.5 .* (1 .- cos.(2π .* (0:nfft-1) ./ (nfft - 1))), 1, nfft, 1)

    windowed_iq = raw_iq .* win
    S = fftshift(fft(windowed_iq, 2), 2)

    spectra_3d = mean(abs.(S) .^ 2, dims=3)
    beam.denoised_spectra = dropdims(spectra_3d, dims=3) ./ num_incoh

    return beam
end

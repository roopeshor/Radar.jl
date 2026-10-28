
"""
    compute_moments(beam)

SIMPLE.COMPUTE_MOMENTS Computes Doppler spectral moments (Woodman method) around peak power.

Locates the maximum spectral power bin for each range bin, traces peak boundaries outward to
adjacent zero-crossings, and calculates zeroth moment (M0, total power), first moment (M1, mean
Doppler shift in Hz), and second central moment (M2, spectral width variance).

# Arguments
- `beam::RadarData`: Radar beam object containing spectra (or denoised_spectra) and timing parameters.

# Returns
- `new_beam::RadarData`: Modified RadarData beam populated with computed M0, M1 (Hz), and M2 profiles.
"""
function compute_moments(beam::RadarData)
    ipp_us = beam.ipp_us
    n_coh = beam.n_coh
    if !isempty(beam.denoised_spectra)
        P_filtered = beam.denoised_spectra
    else
        P_filtered = beam.spectra
    end

    height_bins, nfft = size(P_filtered)

    M0 = zeros(Float64, height_bins)
    M1 = zeros(Float64, height_bins)
    M2 = zeros(Float64, height_bins)
    SNR = zeros(Float64, height_bins)

    for i in 1:height_bins
        P = P_filtered[i, :]
        _, l = findmax(P)

        # Expand peak bounds to local zero-crossings
        minI = l
        maxI = l
        while minI > 1 && P[minI - 1] > 0
            minI -= 1
        end

        while maxI < nfft && P[maxI + 1] > 0
            maxI += 1
        end

        signal_block = P[minI:maxI]
        indices = (minI:maxI) .- 1

        M0[i] = sum(signal_block)
        freq = (indices .- (nfft / 2)) ./ (ipp_us * n_coh * nfft * 1e-6)

        M1[i] = sum(signal_block .* freq) / M0[i]
        M2[i] = sum(((freq .- M1[i]) .^ 2) .* signal_block) / M0[i]
        SNR[i] = sum(((freq .- M1[i]) .^ 3) .* signal_block) / M0[i]
    end

    beam.M0 = M0
    beam.M1 = M1
    beam.M2 = M2
    beam.SNR = SNR

    return beam
end

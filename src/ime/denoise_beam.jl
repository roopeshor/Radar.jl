"""
    denoise_beam(beam)

# Arguments
- `beam::RadarData`: RadarData beam object.

# Returns
- `beam::RadarData`
"""
function denoise_beam(beam::RadarData)
    # 1. Removes ground clutter (DC bin) by interpolation
    # 2. Applies a 5-point moving average filter
    # 3. Estimates and subtracts Hildebrand-Sekhon noise

    Praw = isempty(beam.denoised_spectra) ? beam.spectra : beam.denoised_spectra
    num_bins, nfft = size(Praw)

    dc_idx = (nfft ÷ 2) + 1

    # 1. Ground Clutter Removal — interpolate DC bin and immediate neighbours
    left_val  = Praw[:, dc_idx - 2]
    right_val = Praw[:, dc_idx + 2]
    avg_val   = (left_val .+ right_val) ./ 2
    Praw      = copy(Praw)
    Praw[:, (dc_idx - 1):(dc_idx + 1)] .= avg_val

    # 2. 5-point moving average (manual — replaces MATLAB movmean / missing package)
    kernel   = fill(1.0 / 5, 5)
    pad      = 2   # half-width
    Psmooth  = similar(Praw)
    for i in 1:num_bins
        row = Praw[i, :]
        for j in 1:nfft
            lo  = max(1, j - pad)
            hi  = min(nfft, j + pad)
            Psmooth[i, j] = mean(row[lo:hi])
        end
    end

    # 3. HS noise estimation and subtraction
    M       = beam.m_sNumOfInCohIntegrations
    noise_level = HS_noise_estimate(Psmooth, M)
    Pclean  = Psmooth .- noise_level
    Pclean[Pclean .< 0] .= 0

    beam.denoised_spectra = Pclean
    beam.noise_level      = noise_level

    return beam
end

"""
    compute_beam_spectra_vec(beam; max_axis = 2)

MCCF.compute_beam_spectra Computes power spectra using Mutual Cross-Correlation Function (MCCF).

Processes complex IQ time series by applying frequency-domain DC clutter masking,
segmenting the signal into overlapping windows, and computing cross-correlations
(mutual convolution with Hermitian reversal) across delayed segments to suppress
uncorrelated noise. A Gaussian smoothing filter is applied before FFT shift and
incoherent integration.

# Arguments
- `beam::RadarData`: Radar beam object containing raw IQ complex time-series cube (BeamData).
- `max_axis::Int=2`: which axis to take maximum

# Returns
- `beam::RadarData`
"""
function compute_beam_spectra_vec(beam::RadarData; max_axis::Int = 2)
    beamData = beam.BeamData
    rangeBins, nfft, nInCoh = size(beamData)
    # Total number of signals across all range bins and incoherent integrations
    M = rangeBins * nInCoh
    seg_len = 512
    num_delays = 4
    step = (nfft - seg_len) ÷ (num_delays - 1)

    # Reshape beamData to (M x nfft) where each row is one time-series signal
    signals = reshape(permutedims(beamData, (1, 3, 2)), M, nfft)

    # Ground Clutter Masking around DC bins
    F_sig = fft(signals, 2)
    noise_val = median(abs.(F_sig), dims=2)
    mask_bins = [1, 2, 3, nfft, nfft - 1]
    rand_phase = exp.(im * 2 * pi * rand(M, length(mask_bins)))
    F_sig[:, mask_bins] = noise_val .* rand_phase
    signal_clean = ifft(F_sig, 2)

    # Mutual Convolution via Frequency-Domain Accumulation
    L_conv = 2 * seg_len - 1
    v_r1 = signal_clean[:, 1:seg_len]

    v_r1_padded = zeros(eltype(v_r1), M, L_conv)
    v_r1_padded[:, 1:seg_len] = v_r1
    V_r1_f = fft(v_r1_padded, 2)

    # Accumulate FFT of delayed Hermitian-reversed segments across delays
    V_ri_sum_f = zeros(eltype(v_r1), M, L_conv)
    v_ri_padded = zeros(eltype(v_r1), M, L_conv)

    for d in 2:num_delays
        start_idx = 1 + (d - 1) * step
        v_ri = signal_clean[:, start_idx:(start_idx + seg_len - 1)]
        v_ri_hr = conj(reverse(v_ri, dims=2))
        fill!(v_ri_padded, 0)
        v_ri_padded[:, 1:seg_len] = v_ri_hr
        V_ri_sum_f .+= fft(v_ri_padded, 2)
    end

    # Single inverse FFT to get the accumulated mutual convolution
    conv_accum = ifft(V_r1_f .* V_ri_sum_f, 2) ./ (num_delays - 1)

    # Center-Truncation / Zero-Padding to NFFT length
    padded_sig = zeros(eltype(conv_accum), M, nfft)
    if L_conv < nfft
        offset = (nfft - L_conv) ÷ 2
        padded_sig[:, (offset + 1):(offset + L_conv)] = conv_accum
    elseif L_conv > nfft
        offset = (L_conv - nfft) ÷ 2
        padded_sig = conv_accum[:, (offset + 1):(offset + nfft)]
    else
        padded_sig .= conv_accum
    end

    # Gaussian smoothing filter (5-point window)
    g_kernel = exp.(-((-2:2) .^ 2) ./ 2.0)
    g_kernel ./= sum(g_kernel)

    smoothed_sig = similar(padded_sig)
    for i in 1:M
        row = padded_sig[i, :]
        padded_row = vcat(row[2:-1:1], row, row[end:-1:end-1])
        conv_row = conv(padded_row, g_kernel)
        smoothed_sig[i, :] = conv_row[3:(end - 2)]
    end

    spec_all = abs.(fftshift(fft(smoothed_sig, 2), 2)) .^ 2

    # Reshape back to (rangeBins x nInCoh x nfft) and average over nInCoh (dimension 2)
    spec_3d = reshape(spec_all, rangeBins, nInCoh, nfft)
    spectra = dropdims(mean(spec_3d, dims=2), dims=2)

    max_vals = maximum(spectra, dims=max_axis)
    max_vals[max_vals .== 0] .= 1.0
    beam.denoised_spectra = spectra ./ max_vals

    return beam
end

"""
    compute_mc_spectra(beam; max_axis = 2, apply_gaussian_filt = true)

MCCF.compute_beam_spectra Computes power spectra using Mutual convolution.

Processes complex IQ time series by applying frequency-domain DC clutter masking,
segmenting the signal into overlapping windows, and computing cross-correlations
(mutual convolution with Hermitian reversal) across delayed segments to suppress
uncorrelated noise. A Gaussian smoothing filter is applied before FFT shift and
incoherent integration.

# Arguments
- `beam::RadarData`: Radar beam object containing raw IQ complex time-series cube (BeamData).
- `max_axis`: which axis to take maximum. 'all' would divide whole array by largest spectral peak.
- `apply_gaussian_filt::Bool=true`

# Returns
- RangeBins x NFFT matrix of normalized power spectra.
"""
function compute_mc_spectra(beam::RadarData; max_axis = 2, apply_gaussian_filt::Bool = true)
    beamData = beam.BeamData
    rangeBins, nfft, nInCoh = size(beamData)
    spectraCube = zeros(Float64, rangeBins, nfft, nInCoh)

    seg_len = 512
    num_delays = 4
    step = (nfft - seg_len) ÷ (num_delays - 1)

    for bin in 1:rangeBins
        for inCoh in 1:nInCoh
            signal = beamData[bin, :, inCoh]

            # Frequency-domain ground clutter masking around DC bins
            F_sig = fft(signal)
            noise_val = median(abs.(F_sig))
            # Mask DC and adjacent bins to remove clutter (e.g., bin 1, 2, 3 and N, N-1)
            mask_bins = [1, 2, 3, nfft, nfft - 1]
            # Replace with noise (with random phase to simulate background noise)
            F_sig[mask_bins] = noise_val .* exp.(im * 2 * pi * rand(length(mask_bins)))
            signal_clean = ifft(F_sig)

            # Mutual convolution across delayed segments with Hermitian reversal
            v_r1 = signal_clean[1:seg_len]

            # Accumulator for the convolved results
            conv_accum = zeros(eltype(signal_clean), 2 * seg_len - 1)

            for d in 2:num_delays
                start_idx = 1 + (d - 1) * step
                v_ri = signal_clean[start_idx:(start_idx + seg_len - 1)]
                v_ri_hr = conj(reverse(v_ri))

                v_ci = conv(v_r1, v_ri_hr)
                conv_accum .+= v_ci
            end
            conv_accum ./= (num_delays - 1)

            # Center-truncate or zero-pad accumulated correlation to NFFT length
            final_len = length(conv_accum)
            padded_sig = zeros(eltype(conv_accum), nfft)
            if final_len < nfft
                offset = (nfft - final_len) ÷ 2
                padded_sig[(offset + 1):(offset + final_len)] = conv_accum
            elseif final_len > nfft
                offset = (final_len - nfft) ÷ 2
                padded_sig .= conv_accum[(offset + 1):(offset + final_len)]
            else
                padded_sig .= conv_accum
            end

            # Gaussian adaptive noise filtering prior to Doppler FFT
            smoothed_sig = padded_sig
            if apply_gaussian_filt
                g_kernel = exp.(-((-2:2) .^ 2) ./ 2.0)
                g_kernel ./= sum(g_kernel)
                padded_row = vcat(padded_sig[2:-1:1], padded_sig, padded_sig[end:-1:end-1])
                smoothed_sig = conv(padded_row, g_kernel)[3:(end - 2)]
            end

            spec = abs.(fftshift(fft(smoothed_sig))) .^ 2
            spectraCube[bin, :, inCoh] = spec
        end
    end

    # Incoherent integration and peak normalization
    spectraCube_mean = mean(spectraCube, dims=3)
    spectra = dropdims(spectraCube_mean, dims=3)

    if max_axis == "all" || max_axis == :all
        max_val = maximum(spectra)
        max_val = max_val == 0 ? 1.0 : max_val
        spectra ./= max_val
    else
        max_vals = maximum(spectra, dims=max_axis)
        max_vals[max_vals .== 0] .= 1.0
        spectra ./= max_vals
    end

    return spectra
end

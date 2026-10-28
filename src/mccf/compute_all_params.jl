"""
    compute_all_params(beam; prev_frame_peaks=Int[], apply_mclms=true, win_coeff=2.0, wind_shear=1.0, velocity_res=0.06, cost_weights=[0.5, 0.5])

# Arguments
- `beam::RadarData`: RadarData object with BeamData cube [RangeBins x NFFT x InCoh]
- `prev_frame_peaks::AbstractVector{<:Integer}=Int[]`: Array of peak indices from the previous frame (optional)
- `apply_mclms::Bool=true`: whether to apply mclms filter
- `win_coeff::Real=2.0`: window coefficient
- `wind_shear::Real=1.0`: max wind shear velocity (e.g., 1.0 m/s)
- `velocity_res::Real=0.06`: velocity resolution (e.g., 0.06 m/s)
- `cost_weights::AbstractVector{<:Real}=[0.5, 0.5]`: cost weights: [velocity continuity weight, peak power weight]

# Returns
- `new_beam::RadarData`
"""
function compute_all_params(beam::RadarData;
    prev_frame_peaks::AbstractVector{<:Integer}=Int[],
    apply_mclms::Bool=true,
    win_coeff::Real=2.0,
    wind_shear::Real=1.0,
    velocity_res::Real=0.06,
    cost_weights::AbstractVector{<:Real}=[0.5, 0.5]
)
    beamData = beam.BeamData
    height_bins, nfft, _ = size(beamData)

    # algorithm parameters
    alpha = win_coeff
    max_shear = wind_shear
    dv = velocity_res

    if haskey(beam.a_params, :alpha) && !isnothing(beam.a_params.alpha)
        alpha = beam.a_params.alpha
    end
    if haskey(beam.a_params, :max_shear) && !isnothing(beam.a_params.max_shear)
        max_shear = beam.a_params.max_shear
    end
    if haskey(beam.a_params, :velocity_res) && !isnothing(beam.a_params.velocity_res)
        dv = beam.a_params.velocity_res
    end

    df = round(Int, max_shear / dv) # wind shear bin offset
    w1 = cost_weights[1]
    w2 = cost_weights[2]
    mclms_len = 10                     # mCLMS filter length

    cost_scores = zeros(Float64, height_bins, nfft)
    mclms_denoised_spectra = zeros(Float64, height_bins, nfft)
    peaks_idx = zeros(Int, height_bins)

    F0 = 0 # spectral peak of previous range bin

    spectra = fft(beam.BeamData, 2)
    spectra = mean(spectra, dims=3)
    spectra = dropdims(spectra, dims=3)
    noise_val = median(abs.(spectra), dims=2)
    height_bins, nfft = size(spectra)

    # Zero-frequency DC clutter mask
    mask_bins = [1, 2, 3, nfft, nfft - 1]
    rand_phase = exp.(im * 2 * pi * rand(height_bins, length(mask_bins)))
    spectra[:, mask_bins] = noise_val .* rand_phase
    CS = ifft(spectra ./ maximum(abs.(spectra)), 2)

    # space time doppler window (inclusive bounds)
    # 'height_bins' pairs of (min, max)
    ST_W_array = zeros(Int, height_bins, 2)

    for r in 1:height_bins
        signal_clean = CS[r, :]
        # -----------------------------------------------------------------
        # setting Space-Time Doppler Window (ST-W)
        # -----------------------------------------------------------------
        if r == 1
            # full spectrum for first bin
            ST_W = 1:nfft
        else
            # space window
            W1 = max(1, round(Int, F0-alpha*df)):min(nfft, round(Int, F0+alpha*df))
            if !isempty(prev_frame_peaks) && length(prev_frame_peaks) >= r
                T0 = prev_frame_peaks[r]
                # time window
                W2 = max(1, round(Int, T0-alpha*df)):min(nfft, round(Int, T0+alpha*df))
                inter_W = intersect(W1, W2)
                if length(inter_W) >= alpha * df
                    ST_W = inter_W
                else
                    ST_W = W1
                end
            else
                ST_W = W1
            end
        end
        ST_W_array[r, :] = [minimum(ST_W), maximum(ST_W)]

        # -----------------------------------------------------------------
        # mutual Convolution & Optimal/Suboptimal Vector Selection
        # -----------------------------------------------------------------

        seg_len = nfft ÷ 2
        num_delays = 8
        step = (nfft - seg_len) ÷ (num_delays - 1)

        v_r1 = signal_clean[1:seg_len]
        V_conv = zeros(eltype(signal_clean), num_delays - 1, seg_len)
        max_powers = zeros(Float64, num_delays - 1)

        for d in 2:num_delays
            start_idx = 1 + (d - 1) * step
            v_ri = signal_clean[start_idx:(start_idx+seg_len-1)]
            v_ri_hr = conj(reverse(v_ri))
            v_ci = conv(v_r1, v_ri_hr)

            # get seg_len points from center
            mid_start = seg_len ÷ 2
            v_ci_center = v_ci[mid_start:(mid_start+seg_len-1)]
            V_conv[d-1, :] = v_ci_center

            spec_d = abs.(fftshift(fft(v_ci_center, nfft))) .^ 2
            # max power inside ST-W
            max_powers[d-1] = maximum(spec_d[ST_W])
        end

        _, best_d_idx = findmax(max_powers)
        # optimal solution
        X_opt = V_conv[best_d_idx, :]

        # suboptimal solution is adjacent to optimal
        if best_d_idx > 1
            X_sopt = V_conv[best_d_idx-1, :]
        else
            X_sopt = V_conv[best_d_idx+1, :]
        end

        # -----------------------------------------------------------------
        # modified Complex LMS (MCLMS) Adaptive Filter
        # -----------------------------------------------------------------
        Y_r = X_opt
        if apply_mclms
            Y_r = mclms_filter(X_opt, X_sopt, mclms_len)
        end

        spec_mclms = abs.(fftshift(fft(Y_r, nfft))) .^ 2
        spec_mclms = spec_mclms ./ maximum(spec_mclms)
        mclms_denoised_spectra[r, :] = spec_mclms

        # -----------------------------------------------------------------
        # candidate Peak Selection & Cost Function
        # -----------------------------------------------------------------
        if r == 1
            # initial peak selection
            avg_conv = mean(V_conv, dims=1) |> vec
            final_spectr = abs.(fftshift(fft(avg_conv, nfft))) .^ 2
            _, selected_peak_idx = findmax(final_spectr)
        else
            ST_W_spec = spec_mclms[ST_W]
            M_rmax = maximum(ST_W_spec)

            local_pks = Float64[]
            local_locs = Int[]
            for i in 2:(length(ST_W_spec)-1)
                if ST_W_spec[i] > ST_W_spec[i-1] && ST_W_spec[i] >= ST_W_spec[i+1]
                    push!(local_pks, ST_W_spec[i])
                    push!(local_locs, i)
                end
            end

            pks = local_pks
            locs_in_stw = local_locs

            # peaks that are atleast 40% of max peak
            valid_mask = pks .>= 0.4 * M_rmax
            cand_pks = pks[valid_mask]
            cand_locs = ST_W[locs_in_stw[valid_mask]]

            if isempty(cand_locs)
                _, max_rel_idx = findmax(ST_W_spec)
                cand_locs = [ST_W[max_rel_idx]]
                cand_pks = [M_rmax]
            end

            # select top 3 peaks
            sort_idx = sortperm(cand_pks, rev=true)
            if length(cand_locs) > 3
                cand_locs = cand_locs[sort_idx[1:3]]
                cand_pks = cand_pks[sort_idx[1:3]]
            else
                cand_locs = cand_locs[sort_idx]
                cand_pks = cand_pks[sort_idx]
            end

            best_cost = -Inf
            selected_peak_idx = cand_locs[1]

            for c in 1:length(cand_locs)
                cand_idx = cand_locs[c]
                cand_power = cand_pks[c]

                # velocity cost
                V_cost = 1.0 - abs((cand_idx - F0) * dv) / (alpha * max_shear)

                # cost function
                cost_val = w1 * V_cost + w2 * cand_power
                cost_scores[r, cand_idx] = cost_val

                if cost_val > best_cost
                    best_cost = cost_val
                    selected_peak_idx = cand_idx
                end
            end
        end

        F0 = selected_peak_idx
        peaks_idx[r] = selected_peak_idx

        # -----------------------------------------------------------------
        # moment Estimation (M0, M1, M2)
        # -----------------------------------------------------------------
        l = selected_peak_idx
        minI = l
        maxI = l
        while minI > 1 && spec_mclms[minI-1] > 0
            minI -= 1
        end
        while maxI < nfft && spec_mclms[maxI+1] > 0
            maxI += 1
        end
    end

    fmax = (nfft / 2) / (nfft * beam.ipp_us * 1e-6 * beam.n_coh)
    beam.M1 = -fmax .+ (peaks_idx .- 1.0) .* (fmax - (-fmax)) ./ (1024.0 - 1.0)
    beam.denoised_spectra = mclms_denoised_spectra
    beam.a_params.mccf_cost_scores = cost_scores
    beam.a_params.mccf_ST_W = ST_W_array
    beam.a_params.peaks_idx = peaks_idx

    return beam
end

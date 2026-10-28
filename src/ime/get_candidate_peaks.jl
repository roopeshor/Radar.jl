"""
    get_candidate_peaks(Pclean, noise_level, vel_axis)

IME_ADAPTIVE_WINDOWS (Step B)
Adaptive Doppler-window and prospective peak selection.

# Arguments
- `Pclean::AbstractMatrix{<:Real}`: Cleaned power spectrum matrix.
- `noise_level::AbstractVector{<:Real}`: Noise level vector.
- `vel_axis::AbstractVector{<:Real}`: Velocity axis vector.

# Returns
- Vector of candidate peaks for each range bin.
"""
function get_candidate_peaks(Pclean::AbstractMatrix{<:Real}, noise_level::AbstractVector{<:Real}, vel_axis::AbstractVector{<:Real})
    numRangeBins, nfft = size(Pclean)
    candidates = Vector{Any}(undef, numRangeBins)

    # The Doppler window width is limited to 20% of the entire Doppler axis
    total_bandwidth = maximum(vel_axis) - minimum(vel_axis)
    half_width = (0.20 * total_bandwidth) / 2

    previous_velocity = NaN

    for r in 1:numRangeBins
        P = Pclean[r, :]
        Sn = noise_level[r]

        # 1. Determine Adaptive Window Center
        if r == 1
            _, peak_idx = findmax(P)
            center_velocity = vel_axis[peak_idx]
        else
            center_velocity = previous_velocity
        end

        # 2. Extract the Search Window
        lower_velocity = center_velocity - half_width
        upper_velocity = center_velocity + half_width

        idx_window = findall(v -> v >= lower_velocity && v <= upper_velocity, vel_axis)

        if isempty(idx_window)
            candidates[r] = NamedTuple{(:idx, :velocity, :power, :snr_db, :width), Tuple{Int, Float64, Float64, Float64, Int}}[]
            continue
        end

        Pworking = copy(P[idx_window])
        selected = NamedTuple{(:idx, :velocity, :power, :snr_db, :width), Tuple{Int, Float64, Float64, Float64, Int}}[]

        # 3. Iteratively select up to 3 prospective peaks
        for p in 1:3
            if isempty(Pworking)
                break
            end

            max_amp, local_idx = findmax(Pworking)
            if max_amp <= 0
                break
            end # No more valid signals

            # Find boundaries (zero crossings) within the window
            left = local_idx
            while left > 1 && Pworking[left - 1] > 0
                left -= 1
            end

            right = local_idx
            while right < length(Pworking) && Pworking[right + 1] > 0
                right += 1
            end

            # Extract the actual peak segment data
            P_peak = Pworking[left:right]
            V_peak = vel_axis[idx_window[left:right]]

            # Calculate True Moments
            total_peak_power = sum(P_peak) # 0th Moment (Total Power)
            mean_velocity = sum(P_peak .* V_peak) / total_peak_power # 1st Moment (Power-weighted velocity)

            # Calculate SNR
            if Sn > 0
                total_noise_power = Sn * nfft
                snr_linear = total_peak_power / total_noise_power
            else
                snr_linear = Inf
            end
            snr_db = 10 * log10(snr_linear)

            # Paper requires SNR > 0 dB to be considered a valid peak
            if snr_db <= 0
                break
            end

            # Save candidate
            newCandidate = (
                idx = idx_window[local_idx], # Index of absolute maximum
                velocity = mean_velocity,    # Power-weighted velocity (highly precise)
                power = total_peak_power,
                snr_db = snr_db,
                width = (right - left + 1)
            )

            push!(selected, newCandidate)

            # Mask out this peak with 0s to search for the next one
            Pworking[left:right] .= 0
        end

        candidates[r] = selected

        # 4. Determine center for next range bin (Trace closest peak)
        if !isempty(selected)
            velocities = [c.velocity for c in selected]
            _, closest = findmin(abs.(velocities .- center_velocity))
            previous_velocity = velocities[closest]
        else
            if r == 1
                previous_velocity = center_velocity
            end
        end
    end

    return candidates
end

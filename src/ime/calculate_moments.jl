"""
    calculate_moments(P, vel_axis, noise_level)

Calculates true Doppler spectral moments for a given peak.

# Arguments
- `P::AbstractVector{<:Real}`: Power spectrum vector.
- `vel_axis::AbstractVector{<:Real}`: Velocity axis vector.
- `noise_level::Real`: Estimated noise level.

# Returns
- A named tuple containing `M0`, `M1`, `M2`, `width`, and `SNR`.
"""
function calculate_moments(P::AbstractVector{<:Real}, vel_axis::AbstractVector{<:Real}, noise_level::Real)
    # Noise subtraction and bounding
    Ps = P .- noise_level
    Ps[Ps .< 0] .= 0

    # Find the peak
    max_pow, peak_idx = findmax(Ps)
    if max_pow <= 0
        return (M0 = 0.0, M1 = NaN, M2 = NaN, width = NaN, SNR = -Inf)
    end

    # Trace boundaries of zero-crossings
    left = peak_idx
    while left > 1 && Ps[left - 1] > 0
        left -= 1
    end

    right = peak_idx
    while right < length(Ps) && Ps[right + 1] > 0
        right += 1
    end

    # Extract the continuous signal segment
    idx = left:right
    Pseg = Ps[idx]
    Vseg = vel_axis[idx]

    M0 = sum(Pseg)
    if M0 <= 0
        return (M0 = 0.0, M1 = NaN, M2 = NaN, width = NaN, SNR = -Inf)
    end

    M1 = sum(Vseg .* Pseg) / M0
    M2 = sum(((Vseg .- M1) .^ 2) .* Pseg) / M0
    width = sqrt(M2)

    total_noise = noise_level * length(Ps) # Total noise across the whole NFFT
    SNR = 10 * log10(M0 / total_noise)

    return (M0 = M0, M1 = M1, M2 = M2, width = width, SNR = SNR)
end

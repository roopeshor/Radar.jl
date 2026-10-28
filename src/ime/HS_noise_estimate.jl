"""
    HS_noise_estimate(spectra, in_coh)

# Arguments
- `spectra::AbstractMatrix{<:Real}`
- `in_coh::Real`: incoherent integrations

# Returns
- Vector of noise levels for each range bin.
"""
function HS_noise_estimate(spectra::AbstractMatrix{<:Real}, in_coh::Real)
    numRangeBins, N = size(spectra)
    noise_level = zeros(Float64, numRangeBins)

    if in_coh < 1
        in_coh = 1
    end

    for r in 1:numRangeBins
        A = sort(spectra[r, :])
        n = 1:N
        Pn = cumsum(A) ./ n
        secondMoment = cumsum(A .^ 2) ./ n
        Qn = secondMoment .- Pn .^ 2
        Qn[Qn .<= 0] .= NaN
        Rn = (Pn .^ 2) ./ (in_coh .* Qn)
        valid = findall(x -> x > 1, Rn)
        if isempty(valid)
            # No statistically valid noise region.
            # Conservative fallback: lowest few spectral points.
            k = max(1, round(Int, 0.1 * N))
        else
            k = maximum(valid)
        end
        noise_level[r] = Pn[k]
    end

    return noise_level
end


"""
    HS_noise_estimate(beam)

SIMPLE.HS_NOISE_ESTIMATE Estimates mean noise floor using Hildebrand-Sekhon chi-squared statistics.

Sorts spectral bins in ascending order for each range bin and computes running mean P_n and
running variance Q_n. Scans the test statistic R_n = P_n^2 / (M * Q_n) to isolate the maximum
noise-only subset (R_n > 1) and returns the corresponding mean noise floor.

# Arguments
- `beam::RadarData`: Radar beam object containing spectra and incoherent integration parameters.

# Returns
- Estimated mean noise level vector (1 x RangeBins).
"""
function HS_noise_estimate(beam::RadarData)
    M = beam.m_sNumOfInCohIntegrations
    height_bins, N = size(beam.spectra)

    A = sort(abs.(beam.spectra) .^ 2, dims=2)

    nm1 = reshape((0:(N - 1)) .+ N, 1, N)
    P_n = cumsum(A, dims=2) ./ nm1
    nm2 = reshape(1:N, 1, N)
    Q_n = (cumsum(A .^ 2, dims=2) ./ nm2) .- (P_n .^ 2)

    R_n = (P_n .^ 2) ./ (Q_n .* M)

    noise_level = zeros(Float64, height_bins)

    for i in 1:height_bins
        arr = reverse(R_n[i, :])

        idx_0based = 0
        for j in 1:length(arr)
            if !isnan(arr[j]) && arr[j] > 1.0
                idx_0based = j - 1
                break
            end
        end

        valid_idx_0based = N - 1 - idx_0based
        noise_level[i] = P_n[i, valid_idx_0based + 1]
    end

    return noise_level
end

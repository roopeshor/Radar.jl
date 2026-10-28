"""
    synthesize_iq_data_vectrzd(; vr, Header, spec_w, SNR=[], add_iq_noise=true, dc_rel_power=0.0, gc_rel_power=0.005, sigma_dc=[], sigma_gc=[], bright_bands=Matrix{Float64}(undef, 0, 2))

Synthesizes complex IQ time-series data using Zrnic spectral simulation.

   Generates Gaussian power spectrum signatures modeled on Doppler velocity, spectral width, and SNR,
   applies complex Gaussian noise shaping (Zrnic method), and computes inverse FFT to yield IQ time series.
Vectorized function

# Arguments
- `vr::AbstractVector{<:Real}`: Radial velocity profile vector (m/s).
- `Header::BeamHeader`: BeamHeader object containing NFFT, RangeBins, IPP, and coherent integration parameters.
- `spec_w::AbstractVector{<:Real}`: Spectral width (m/s).
- `SNR::AbstractVector{<:Real}=Float64[]`: Signal-to-noise ratio array across range bins.
- `add_iq_noise::Bool=true`: Flag to enable complex IQ Gaussian noise synthesis.
- `dc_rel_power::Real=0.0`: additional dc offset. this fraction of peak value of spectrum is added to signal in time domain
- `gc_rel_power::Real=0.005`: ground clutter power. this fraction of peak value of spectrum is added to signal in time domain
- `sigma_dc::AbstractVector{<:Real}=Float64[]`: width of DC component across range bins.
- `sigma_gc::AbstractVector{<:Real}=Float64[]`: width of ground clutter across range bins.
- `bright_bands::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 2)`: list of height bin and bright band widths

# Returns
- `BeamData::Array{ComplexF64, 3}`: Synthesized complex IQ time-series data cube (RangeBins x NFFT x InCoh).
"""
function synthesize_iq_data_vectrzd(;
    vr::AbstractVector{<:Real},
    Header::BeamHeader,
    spec_w::AbstractVector{<:Real},
    SNR::AbstractVector{<:Real}=Float64[],
    add_iq_noise::Bool=true,
    dc_rel_power::Real=0.0,
    gc_rel_power::Real=0.005,
    sigma_dc::AbstractVector{<:Real}=Float64[],
    sigma_gc::AbstractVector{<:Real}=Float64[],
    bright_bands::AbstractMatrix{<:Real}=Matrix{Float64}(undef, 0, 2)
)
    NFFT = Header.m_sNFFT
    nRangeBins = Header.m_sNumOfRangeBins
    ipp_us = Header.m_fIntrPulsePeriod_us
    cohIn = Header.m_sNumOfCohIntegrations

    P_noise = 1000.0
    wavelen = 3.0e8 / 205e6
    DBS_V = wavelen / 2.0
    fd = .-vr ./ DBS_V
    PRF_eff = 1.0 / (ipp_us * 1e-6 * cohIn)
    sigma_f = 2.0 .* spec_w ./ wavelen

    SNR_vec = isempty(SNR) ? (10.0 .^ (range(40, -40, length=nRangeBins) ./ 10.0)) : vec(SNR)
    sigma_dc_vec = isempty(sigma_dc) ? fill(1.0 / 20.0, nRangeBins) : vec(sigma_dc)
    sigma_gc_vec = isempty(sigma_gc) ? (1.0 ./ ((1:nRangeBins) ./ 10.0) ./ 5.0) : vec(sigma_gc)

    f = range(-PRF_eff / 2.0, PRF_eff / 2.0 - PRF_eff / NFFT, length=NFFT)
    bright_band_ranges = isempty(bright_bands) ? Int[] : round.(Int, bright_bands[:, 1])

    fd_col = vec(fd)
    sigma_f_col = vec(sigma_f)
    SNR_col = vec(SNR_vec)
    sigma_dc_col = vec(sigma_dc_vec)
    sigma_gc_col = vec(sigma_gc_vec)
    f_row = reshape(f, 1, NFFT)

    P_sig = P_noise .* SNR_col
    S_sig = exp.(-((f_row .- fd_col) .^ 2) ./ (2.0 .* sigma_f_col .^ 2))
    sum_S_sig = sum(S_sig, dims=2)

    scale = zeros(Float64, nRangeBins)
    valid = vec(sum_S_sig) .> 0.0
    scale[valid] = (P_sig[valid] ./ vec(sum_S_sig)[valid]) ./ (sqrt(2.0 * pi) * sigma_f_col[valid])
    S_sig = S_sig .* scale

    S = S_sig .+ (P_noise / NFFT)

    # DC Component
    sum_S = sum(S, dims=2)
    S = S .+ (dc_rel_power .* sum_S) .* exp.(-(f_row .^ 2) ./ (2.0 .* sigma_dc_col .^ 2))

    # Ground Clutter
    sum_S = sum(S, dims=2)
    S = S .+ (gc_rel_power .* sum_S) .* exp.(-(f_row .^ 2) ./ (2.0 .* sigma_gc_col .^ 2))

    # Bright Bands
    if !isempty(bright_bands)
        num_bb = size(bright_bands, 1)
        for k in 1:num_bb
            z_idx = bright_band_ranges[k]
            if z_idx >= 1 && z_idx <= nRangeBins
                rescaled_val = -3.5 + (k - 1) * (2.0 - (-3.5)) / max(1, (nRangeBins - 1))
                ampl = sum(S[z_idx, :]) * 10.0 ^ rescaled_val
                S[z_idx, :] .+= ampl .* exp.(-(f_row .^ 2) / (2.0 * bright_bands[k, 2]^2))
            end
        end
    end

    S_sqrt = sqrt.(S ./ 2.0)
    Z = if add_iq_noise
        X = randn(nRangeBins, NFFT)
        Y = randn(nRangeBins, NFFT)
        S_sqrt .* (X .+ im .* Y)
    else
        S_sqrt
    end

    BeamData = round.(ifft(ifftshift(Z, 2), 2) .* sqrt(NFFT))
    return reshape(BeamData, Int(nRangeBins), Int(NFFT), 1)  # (RangeBins × NFFT × 1 InCoh)
end

"""
    compute_spectrum_broadening(z, u, v, w, v_radial; sigma_turb_base=0.2, sigma_chi=0.2)

# Arguments
- `z::AbstractVector{<:Real}`
- `u::AbstractVector{<:Real}`
- `v::AbstractVector{<:Real}`
- `w::AbstractVector{<:Real}`
- `v_radial::AbstractVector{<:Real}`
- `sigma_turb_base::Real=0.2`: base turbulence
- `sigma_chi::Real=0.2`: Intermittency severity parameter

# Returns
- `sigma_v::Vector{Float64}`
"""
function compute_spectrum_broadening(z::AbstractVector{<:Real},
    u::AbstractVector{<:Real},
    v::AbstractVector{<:Real},
    w::AbstractVector{<:Real},
    v_radial::AbstractVector{<:Real};
    sigma_turb_base::Real=0.2,
    sigma_chi::Real=0.2
)::Vector{Float64}
    N_gates = length(z)
    dz_m = (z[2] - z[1]) * 1000 # Range gate resolution in meters

    # cross-beam velocity magnitude
    v_cross = sqrt.(u .^ 2 .+ v .^ 2 .+ w .^ 2 .- v_radial .^ 2)

    # Vertical shear broadening across gate width
    dv_dz = abs.(calc_gradient(v_radial, dz_m))
    sigma_shear = (dv_dz .* dz_m) / sqrt(12)

    # Antenna beam broadening (finite one-way beamwidth)
    theta_3db = deg2rad(3.3)
    sigma_beam = (theta_3db / (2 * sqrt(2 * log(2)))) * abs.(v_cross)

    # turbulence
    turb_base_vec = fill(Float64(sigma_turb_base), N_gates) # [m/s]

    # Unmodulated spectral width
    sigma_v_base = sqrt.(sigma_shear .^ 2 .+ sigma_beam .^ 2 .+ turb_base_vec .^ 2)

    # Spatially Correlated Log-Normal Stochastic Process
    mu_chi = -0.5 * sigma_chi^2   # Mean-preserving offset

    raw_noise = randn(N_gates)

    # Apply Gaussian spatial smoothing filter (correlation length ~ 450 m)
    corr_length_gates = round(Int, 0.45 / (z[2] - z[1]))
    window_size = max(3, corr_length_gates)
    eta_corr = smooth_gaussian(raw_noise, window_size)
    eta_corr ./= std(eta_corr) # Renormalize variance to 1

    # Compute log-normal multiplier
    chi = exp.(mu_chi .+ sigma_chi .* eta_corr)

    # final spectrum broadening Width
    sigma_v = sigma_v_base .* chi

    return sigma_v
end

function calc_gradient(y::AbstractVector{<:Real}, h::Real)
    n = length(y)
    g = similar(y, Float64)
    if n == 1
        g[1] = 0.0
        return g
    end
    g[1] = (y[2] - y[1]) / h
    for i in 2:(n-1)
        g[i] = (y[i+1] - y[i-1]) / (2 * h)
    end
    g[n] = (y[n] - y[n-1]) / h
    return g
end

function smooth_gaussian(x::AbstractVector{<:Real}, window_len::Integer)
    half_w = window_len ÷ 2
    sigma = window_len / 4.0
    idx = (-half_w):half_w
    kernel = exp.(-(idx .^ 2) ./ (2 * sigma^2))
    kernel ./= sum(kernel)

    n = length(x)
    padded = vcat(fill(x[1], half_w), x, fill(x[end], half_w))
    smoothed = conv(padded, kernel)
    start_idx = half_w + 1
    return smoothed[start_idx:(start_idx+n-1)]
end

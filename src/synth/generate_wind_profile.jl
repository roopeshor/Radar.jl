"""
# Arguments
- `nRangeBins::Real=173`: [ ] number of range bins
- `start_height::Real=6.030`: [km] max height for which wind profile has to be generated
- `end_height::Real=32.0`: [km] max height for which wind profile has to be generated
- `bg_z0_m::Real=0.05`: [m] Surface roughness length for terrain.
- `bg_ustar::Real=0.3`: [m/s] Friction velocity representing surface stress.
- `bg_kappa::Real=0.4`: [ ] von Kármán constant.
- `bg_blh::Real=1.5`: [km] boundary layer height
- `sj_w::Real=12.0`: [m/s] (zonal) onshore wind speeds of the SW monsoon jet
- `sj_uv::Real=14.0`: [m/s] (meridional) onshore wind speeds of the SW monsoon jet
- `sj_center::Real=1.5`: [km] Altitude of the jet core center.
- `sj_s::Real=0.8`: [km] Gaussian standard deviation controlling vertical thickness of the jet.
- `tej_u::Real=20.0`: [m/s] (zonal) Strong easterly (westward) jet core speed near the tropopause.
- `tej_v::Real=3.0`: [m/s] (zonal) Strong easterly (westward) jet core speed near the tropopause.
- `tej_alt::Real=15.0`: [km] TEJ core altitude
- `tej_vthick::Real=2.5`: [km] TEJ core vertical layer thickness.
- `qbo_lambda::Real=20.0`: [km] Vertical wavelength of alternating wind regimes in the stratosphere.
- `qbo_u::Real=8.0`: [m/s] Peak zonal amplitude of QBO phase shifts.
- `qbo_v::Real=3.0`: [m/s] Peak meridional amplitude of QBO phase shifts.
- `qbo_alt::Real=17.0`: [km] Altitude of QBO region
- `ou_alpha_slope::Real=atan(3/20.0)`: [ ] Western Ghats ridge terrain slope (1.5 km height rise over 30 km distance).
- `ou_H_lift::Real=2.0`: [km] Vertical decay scale
- `gw_nmodes::Real=20.0`: [ ] Number of continuous gravity wave spectrum modes
- `gw_amp::Real=2.5`: [m/s] Scaling factor for gravity wave horizontal velocity
- `gw_max_lambda::Real=4.0`: [km] Upper limit for gravity wave vertical scale
- `gw_damp_scale::Real=12.0`: [km] damping scale to control upper amplitude
- `gw_cf::Real=0.05`: [ ] Coupling factor deriving vertical wave speed from horizontal wave speed
- `lw_Hs::Real=7.0`: [km] Density scale height for wave amplification (e^(z / 2H_s))
- `jitter_amp::Real=0.4`: [m/s] Amplitude of low-frequency spatial jet jitter
- `plume_count::Real=4.0`: [ ] Number of random convective updrafts/downdrafts
- `plume_max_alt::Real=6.0`: [km] Maximum altitude ceiling for thermal plumes
- `kt_vaf::Real=0.25`: [ ] Vertical anisotropy scaling factor
- `kt_bf::Real=1.0`: [ ] Baseline factor for turbulence scaling
- `kt_maxT::Real=2.5`: [ ] Maximum turbulence boost under instability
- `kt_Ric::Real=0.25`: [ ] Critical Richardson Number Threshold
- `kt_tlpeak::Real=3.0`: [km] Turbulence layer peak altitude
- `kt_vthick::Real=2.5`: [km] Vertical thickness scale of turbulent layer
- `lw_N_bv::Real=0.012`: [rad/s] Brunt-Väisälä buoyancy frequency
- `turb_interm_sigma::Real=0.5`: [ ] Intermittency parameter for log-normal weighting

# Returns
- Tuple of `(z, u, v, w)`
"""
function generate_wind_profile(;
    nRangeBins::Real=173,
    start_height::Real=6.030,
    end_height::Real=32.0,
    bg_z0_m::Real=0.05,
    bg_ustar::Real=0.3,
    bg_kappa::Real=0.4,
    bg_blh::Real=1.5,
    sj_w::Real=12.0,
    sj_uv::Real=14.0,
    sj_center::Real=1.5,
    sj_s::Real=0.8,
    tej_u::Real=20.0,
    tej_v::Real=3.0,
    tej_alt::Real=15.0,
    tej_vthick::Real=2.5,
    qbo_lambda::Real=20.0,
    qbo_u::Real=8.0,
    qbo_v::Real=3.0,
    qbo_alt::Real=17.0,
    ou_alpha_slope::Real=atan(3 / 20.0),
    ou_H_lift::Real=2.0,
    gw_nmodes::Real=20.0,
    gw_amp::Real=2.5,
    gw_max_lambda::Real=4.0,
    gw_damp_scale::Real=12.0,
    gw_cf::Real=0.05,
    lw_Hs::Real=7.0,
    jitter_amp::Real=0.4,
    plume_count::Real=4.0,
    plume_max_alt::Real=6.0,
    kt_vaf::Real=0.25,
    kt_bf::Real=1.0,
    kt_maxT::Real=2.5,
    kt_Ric::Real=0.25,
    kt_tlpeak::Real=3.0,
    kt_vthick::Real=2.5,
    lw_N_bv::Real=0.012,
    turb_interm_sigma::Real=0.5
)
    # Height grid configuration (km)
    h = range(start_height, end_height, length=Int(nRangeBins)) |> collect
    N = length(h)

    # Synoptic Background Profile
    u_bg = synoptic_bg(h; bg_z0_m=bg_z0_m, bg_ustar=bg_ustar, bg_kappa=bg_kappa, bg_blh=bg_blh)

    # Low-Level Somali Jet
    u_llj, v_llj = somali_jet(h; sj_w=sj_w, sj_uv=sj_uv, sj_center=sj_center, sj_s=sj_s)

    # Tropical Easterly Jet
    u_tej, v_tej = tropical_easterly_jet(h; tej_u=tej_u, tej_v=tej_v, tej_alt=tej_alt, tej_vthick=tej_vthick)

    # Stratospheric QBO Regime
    u_qbo, v_qbo = stratospheric_qbo(h; qbo_lambda=qbo_lambda, qbo_u=qbo_u, qbo_v=qbo_v, qbo_alt=qbo_alt)

    # Horizontal components with Spatial Jitter
    u_mean = u_llj .+ u_bg .+ u_tej .+ u_qbo
    v_mean = v_llj .+ v_tej .+ v_qbo

    jitter_u = jitter_amp .* smooth_gaussian(randn(N), 30)'
    jitter_v = jitter_amp .* smooth_gaussian(randn(N), 30)'
    u_mean = u_mean .+ vec(jitter_u)
    v_mean = v_mean .+ vec(jitter_v)

    # Orographic Updraft
    w_oro = orographic_updraft(h, u_mean; ou_alpha_slope=ou_alpha_slope, ou_H_lift=ou_H_lift)

    # Continuous Gravity Wave Spectrum
    u_gw, v_gw, w_gw = gravity_waves(h; gw_nmodes=gw_nmodes, gw_amp=gw_amp, gw_max_lambda=gw_max_lambda, gw_damp_scale=gw_damp_scale, gw_cf=gw_cf, lw_Hs=lw_Hs)

    # Random Convective Plumes in Lower Troposphere
    w_plumes = convective_plumes(h; plume_count=Int(plume_count), plume_max_alt=plume_max_alt)

    # Stable Kolmogorov Turbulence with Log-Normal Intermittency
    u_turb, v_turb, w_turb = kolmogorov_turbulence(h, u_mean, v_mean;
        kt_vaf=kt_vaf, kt_bf=kt_bf, kt_maxT=kt_maxT, kt_Ric=kt_Ric, kt_tlpeak=kt_tlpeak, kt_vthick=kt_vthick, turb_interm_sigma=turb_interm_sigma, lw_N_bv=lw_N_bv)

    # 3D Wind Vector Assembly
    z = h
    u = vec(u_mean .+ u_gw .+ u_turb)
    v = vec(v_mean .+ v_gw .+ v_turb)
    w = vec(w_oro .+ w_gw .+ w_plumes .+ w_turb)

    return z, u, v, w
end

"""
    convective_plumes(h; plume_count=4, plume_max_alt=6.0)

# Arguments
- `h::AbstractVector{<:Real}`: height bins
- `plume_count::Integer=4`: Number of random convective updrafts/downdrafts
- `plume_max_alt::Real=6.0`: Maximum altitude ceiling for thermal plumes

# Returns
- `w_plumes::Vector{Float64}`
"""
function convective_plumes(h::AbstractVector{<:Real}; plume_count::Integer=4, plume_max_alt::Real=6.0)::Vector{Float64}
    # Random Convective Plumes in Lower Troposphere
    z = vec(h)
    N = length(z)
    w_plumes_col = zeros(Float64, N)
    for p in 1:plume_count
        z_center = 0.8 + rand() * (plume_max_alt - 0.8)
        w_peak = (rand() - 0.3) * 2.5
        w_width = 0.2 + rand() * 0.4
        w_plumes_col .+= w_peak .* exp.(-((z .- z_center) ./ w_width) .^ 2)
    end
    w_plumes = w_plumes_col'
    return vec(w_plumes)
end

"""
    gravity_waves(h; gw_nmodes=20, gw_amp=2.5, gw_max_lambda=4.0, gw_damp_scale=12.0, gw_cf=0.05, lw_Hs=7.0)

# Arguments
- `h::AbstractVector{<:Real}`: height bins
- `gw_nmodes::Integer=20`: Number of continuous gravity wave spectrum modes
- `gw_amp::Real=2.5`: Scaling factor for gravity wave horizontal velocity
- `gw_max_lambda::Real=4.0`: Upper limit for gravity wave vertical scale
- `gw_damp_scale::Real=12.0`: Tightened damping scale (was 18.0 km) to control upper amplitude
- `gw_cf::Real=0.05`: Coupling factor deriving vertical wave speed from horizontal wave speed
- `lw_Hs::Real=7.0`: Density scale height for wave amplification (e^(z / 2H_s))

# Returns
- Tuple of `(u_gw, v_gw, w_gw)`
"""
function gravity_waves(h::AbstractVector{<:Real};
    gw_nmodes::Integer=20,
    gw_amp::Real=2.5,
    gw_max_lambda::Real=4.0,
    gw_damp_scale::Real=12.0,
    gw_cf::Real=0.05,
    lw_Hs::Real=7.0
)
    # Continuous Gravity Wave Spectrum
    h_col = vec(h)
    N = length(h_col)

    m_wavenumbers = range(2 * pi / gw_max_lambda, 2 * pi / 0.4, length=gw_nmodes) |> collect
    phases_u = rand(gw_nmodes) .* 2 .* pi
    phases_v = rand(gw_nmodes) .* 2 .* pi

    E_m = (m_wavenumbers) .^ (-1.5)
    E_m ./= sum(E_m)

    u_wave_spec = zeros(Float64, N)
    v_wave_spec = zeros(Float64, N)
    for i in 1:gw_nmodes
        u_wave_spec .+= sqrt(E_m[i]) .* cos.(m_wavenumbers[i] .* h_col .+ phases_u[i])
        v_wave_spec .+= sqrt(E_m[i]) .* sin.(m_wavenumbers[i] .* h_col .+ phases_v[i])
    end

    amp_growth = exp.(h_col ./ (2 * lw_Hs)) .* exp.(-h_col ./ gw_damp_scale)
    u_gw = (gw_amp .* (u_wave_spec ./ std(u_wave_spec)) .* amp_growth)'
    v_gw = (gw_amp .* (v_wave_spec ./ std(v_wave_spec)) .* amp_growth)'
    w_gw = -gw_cf .* u_gw

    return u_gw, v_gw, w_gw
end

"""
    kolmogorov_turbulence(h, u_mean, v_mean; kt_vaf=0.25, kt_bf=1.0, kt_maxT=2.5, kt_Ric=0.25, kt_tlpeak=3.0, kt_vthick=2.5, turb_interm_sigma=0.5, lw_N_bv=0.012)

Stable Kolmogorov Turbulence with Log-Normal Intermittency

# Arguments
- `h::AbstractVector{<:Real}`: height bins
- `u_mean::AbstractVector{<:Real}`: mean zonal wind profile
- `v_mean::AbstractVector{<:Real}`: mean meridional wind profile
- `kt_vaf::Real=0.25`: Vertical anisotropy scaling factor
- `kt_bf::Real=1.0`: Baseline factor for turbulence scaling
- `kt_maxT::Real=2.5`: Maximum turbulence boost under instability
- `kt_Ric::Real=0.25`: Critical Richardson Number Threshold
- `kt_tlpeak::Real=3.0`: Turbulence layer peak altitude
- `kt_vthick::Real=2.5`: Vertical thickness scale of turbulent layer
- `turb_interm_sigma::Real=0.5`: Intermittency parameter for log-normal weighting
- `lw_N_bv::Real=0.012`: Brunt-Väisälä buoyancy frequency

# Returns
- Tuple of `(u_turb, v_turb, w_turb)`
"""
function kolmogorov_turbulence(h::AbstractVector{<:Real},
    u_mean::AbstractVector{<:Real},
    v_mean::AbstractVector{<:Real};
    kt_vaf::Real=0.25,
    kt_bf::Real=1.0,
    kt_maxT::Real=2.5,
    kt_Ric::Real=0.25,
    kt_tlpeak::Real=3.0,
    kt_vthick::Real=2.5,
    turb_interm_sigma::Real=0.5,
    lw_N_bv::Real=0.012
)
    z = vec(h)
    N = length(z)
    dz = z[2] - z[1]

    u_m = vec(u_mean)
    v_m = vec(v_mean)

    du_dz = calc_gradient(u_m, dz * 1000.0)
    dv_dz = calc_gradient(v_m, dz * 1000.0)
    shear_sq = du_dz .^ 2 .+ dv_dz .^ 2 .+ 1e-8

    Ri = (lw_N_bv^2) ./ shear_sq
    turb_boost = kt_bf .+ kt_maxT .* (Ri .< kt_Ric) .* exp.(-((z .- kt_tlpeak) ./ kt_vthick) .^ 2)

    k_vec = vcat(0:floor(Int, N/2), (-floor(Int, (N-1)/2)):-1) ./ (N * dz)
    k_min = 1.0 / (N * dz)
    k_abs = max.(abs.(k_vec), k_min)

    H_turb = k_abs .^ (-5 / 6)
    H_turb[1] = 0.0

    raw_u = real(ifft(fft(randn(N)) .* H_turb))
    raw_v = real(ifft(fft(randn(N)) .* H_turb))
    raw_w = real(ifft(fft(randn(N)) .* H_turb))

    norm_u = (raw_u .- mean(raw_u)) ./ std(raw_u)
    norm_v = (raw_v .- mean(raw_v)) ./ std(raw_v)
    norm_w = (raw_w .- mean(raw_w)) ./ std(raw_w)

    mu_interm = -0.5 * turb_interm_sigma^2
    interm_weight = exp.(mu_interm .+ turb_interm_sigma .* norm_u)

    max_sqrt_shear = maximum(sqrt.(shear_sq))
    sigma_base = 0.15 .+ 0.4 .* (sqrt.(shear_sq) ./ max_sqrt_shear)
    sigma_turb = sigma_base .* turb_boost .* interm_weight

    u_turb = (norm_u .* sigma_turb)'
    v_turb = (norm_v .* sigma_turb)'
    w_turb = (norm_w .* (sigma_turb .* kt_vaf))'

    return vec(u_turb), vec(v_turb), vec(w_turb)
end

"""
    orographic_updraft(h, u_mean; ou_alpha_slope=atan(3/20.0), ou_H_lift=2.0)

Orographic Updraft

# Arguments
- `h::AbstractVector{<:Real}`: height bins
- `u_mean::AbstractVector{<:Real}`: mean zonal wind profile
- `ou_alpha_slope::Real=atan(3/20.0)`: [ ] Western Ghats ridge terrain slope (1.5 km height rise over 30 km distance).
- `ou_H_lift::Real=2.0`: [km] Vertical decay scale

# Returns
- `w_oro::Vector{Float64}`
"""
function orographic_updraft(h::AbstractVector{<:Real},
    u_mean::AbstractVector{<:Real};
    ou_alpha_slope::Real=atan(3 / 20.0),
    ou_H_lift::Real=2.0
)
    # Orographic Updraft
    w_oro = max.(0.0, abs.(u_mean)) .* tan(ou_alpha_slope) .* exp.(-h ./ ou_H_lift)
    return vec(w_oro)
end

"""
    somali_jet(h; sj_w=12.0, sj_uv=14.0, sj_center=1.5, sj_s=0.8)

Low-Level Somali Jet (~1.5 km peak) - Onshore SW Monsoon flow

# Arguments
- `h::AbstractVector{<:Real}`: height bins
- `sj_w::Real=12.0`: [m/s] (zonal) onshore wind speeds of the SW monsoon jet
- `sj_uv::Real=14.0`: [m/s] (meridional) onshore wind speeds of the SW monsoon jet
- `sj_center::Real=1.5`: [km] Altitude of the jet core center.
- `sj_s::Real=0.8`: [km] Gaussian standard deviation controlling vertical thickness of the jet.

# Returns
- Tuple of `(u_llj, v_llj)`
"""
function somali_jet(h::AbstractVector{<:Real};
    sj_w::Real=12.0,
    sj_uv::Real=14.0,
    sj_center::Real=1.5,
    sj_s::Real=0.8
)
    # Low-Level Somali Jet (~1.5 km peak) - Onshore SW Monsoon flow
    h_vec = vec(h)
    u_llj = sj_w .* exp.(-((h_vec .- sj_center) ./ sj_s) .^ 2)
    v_llj = sj_uv .* exp.(-((h_vec .- sj_center) ./ sj_s) .^ 2)
    return vec(u_llj), vec(v_llj)
end

"""
    stratospheric_qbo(h; qbo_lambda=20.0, qbo_u=8.0, qbo_v=3.0, qbo_alt=17.0)

Stratospheric QBO Regime (17 - 32 km)

# Arguments
- `h::AbstractVector{<:Real}`: height bins
- `qbo_lambda::Real=20.0`: [km] Vertical wavelength of alternating wind regimes in the stratosphere.
- `qbo_u::Real=8.0`: [m/s] Peak zonal amplitude of QBO phase shifts.
- `qbo_v::Real=3.0`: [m/s] Peak meridional amplitude of QBO phase shifts.
- `qbo_alt::Real=17.0`: [km] Altitude of QBO region

# Returns
- Tuple of `(u_qbo, v_qbo)`
"""
function stratospheric_qbo(h::AbstractVector{<:Real};
    qbo_lambda::Real=20.0,
    qbo_u::Real=8.0,
    qbo_v::Real=3.0,
    qbo_alt::Real=17.0
)
    # Stratospheric QBO Regime (17 - 32 km)
    h_vec = vec(h)
    qbo_envelope = 0.5 .* (1.0 .+ tanh.((h_vec .- qbo_alt) ./ 2.0))
    u_qbo = qbo_u .* sin.(2 * pi .* (h_vec .- qbo_alt) ./ qbo_lambda) .* qbo_envelope
    v_qbo = qbo_v .* cos.(2 * pi .* (h_vec .- qbo_alt) ./ qbo_lambda) .* qbo_envelope

    return vec(u_qbo), vec(v_qbo)
end

"""
    synoptic_bg(h; bg_z0_m=0.05, bg_ustar=0.3, bg_kappa=0.4, bg_blh=1.5)

Synoptic Background Profile
Smooth logarithmic boundary layer fading by 2 km

# Arguments
- `h::AbstractVector{<:Real}`: height bins
- `bg_z0_m::Real=0.05`: [m] Surface roughness length for terrain.
- `bg_ustar::Real=0.3`: [m/s] Friction velocity representing surface stress.
- `bg_kappa::Real=0.4`: [ ] von Kármán constant.
- `bg_blh::Real=1.5`: [km] boundary layer height

# Returns
- `u_bl::Vector{Float64}`
"""
function synoptic_bg(h::AbstractVector{<:Real};
    bg_z0_m::Real=0.05,
    bg_ustar::Real=0.3,
    bg_kappa::Real=0.4,
    bg_blh::Real=1.5
)
    # Synoptic Background Profile
    # Smooth logarithmic boundary layer fading by 2 km
    z_m = h .* 1000.0
    u_bl = (bg_ustar / bg_kappa) .* log.((z_m .+ bg_z0_m) ./ bg_z0_m) .* exp.(-h ./ bg_blh)

    return vec(u_bl)
end

"""
    tropical_easterly_jet(h; tej_u=10.0, tej_v=3.0, tej_alt=15.0, tej_vthick=2.5)

Tropical Easterly Jet (~15 km peak)

# Arguments
- `h::AbstractVector{<:Real}`: height bins
- `tej_u::Real=10.0`: [m/s] (zonal) Strong easterly (westward) jet core speed near the tropopause.
- `tej_v::Real=3.0`: [m/s] (zonal) Strong easterly (westward) jet core speed near the tropopause.
- `tej_alt::Real=15.0`: [km] TEJ core altitude
- `tej_vthick::Real=2.5`: [km] TEJ core vertical layer thickness.

# Returns
- Tuple of `(u_tej, v_tej)`
"""
function tropical_easterly_jet(h::AbstractVector{<:Real};
    tej_u::Real=10.0,
    tej_v::Real=3.0,
    tej_alt::Real=15.0,
    tej_vthick::Real=2.5
)
    # Tropical Easterly Jet (~15 km peak)
    h_vec = vec(h)
    u_tej = tej_u .* exp.(-((h_vec .- tej_alt) ./ tej_vthick) .^ 2)
    v_tej = tej_v .* exp.(-((h_vec .- tej_alt) ./ tej_vthick) .^ 2)

    return vec(u_tej), vec(v_tej)
end

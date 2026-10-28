"""
synth.create_synthetic_observation Synthesizes an Observation object with ground-truth wind profiles.

Generates 5-beam IQ data cubes via Zrnic spectral simulation, sets reference moments and wind vector
profiles (u, v, w), and optionally writes out a .raw binary data file.

# Arguments
- `u::AbstractVector{<:Real}`: Zonal wind profile vector (m/s).
- `v::AbstractVector{<:Real}`: Meridional wind profile vector (m/s).
- `w::AbstractVector{<:Real}`: Vertical wind profile vector (m/s).
- `headerFields::AbstractVector`: Custom header property overrides.
- `filepath::AbstractString=""`: Destination filepath to save synthetic .raw binary file.
- `sigma_v::Any=""`: spectrum broadening
- `SNR::AbstractVector{<:Real}=Float64[]`: SNR array across range bins.
- `sigma_dc::AbstractVector{<:Real}=Float64[]`: width of DC component across range bins.
- `add_iq_noise::Bool=true`: Flag to add complex Gaussian noise to IQ series.
- `dc_rel_power::Real=0.0`: additional dc offset. this fraction of peak value of spectrum is added to signal in time domain
- `gc_rel_power::Real=0.005`: ground clutter power. this fraction of peak value of spectrum is added to signal in time domain
- `bright_bands::Any=NamedTuple()`: list of height bin and bright band widths per beam
- `sigma_turb_base::Real=0.0`
- `sigma_chi::Real=0.1`

# Returns
- `obj::Observation`: Synthesized Observation object populated with 5 beams and reference profiles.
"""
function create_synthetic_observation(u::AbstractVector{<:Real}, v::AbstractVector{<:Real}, w::AbstractVector{<:Real};
    headerFields::AbstractVector,
    filepath::AbstractString="",
    sigma_v::Any="",
    SNR::AbstractVector{<:Real}=Float64[],
    sigma_dc::AbstractVector{<:Real}=Float64[],
    add_iq_noise::Bool=true,
    dc_rel_power::Real=0.0,
    gc_rel_power::Real=0.005,
    bright_bands::Any=NamedTuple(),
    sigma_turb_base::Real=0.0,
    sigma_chi::Real=0.1
)

    Headers = Vector{RadarData}(undef, 5)
    H = headerFields[1]
    range_res = 3.0e8 * (H.m_fBaudLength_us * 1e-6) / 2.0
    z = H.m_fWindow1StartHeight .+ (0:(H.m_sNumOfRangeBins-1)) .* range_res

    DBS_V = 3.0e8 / 205e6 / 2.0

    for k in 1:5
        H = headerFields[k]
        Headers[k] = RadarData(H)
        vr = compute_radial_velocity(u, v, w, H.m_fAzimuth, H.m_fOffZenith)
        sv = sigma_v
        if sv isa AbstractString && sv == ""
            sv = compute_spectrum_broadening(z, u, v, w, vr, sigma_turb_base, sigma_chi)
        end
        Headers[k].ref_M1 = .-vr ./ DBS_V
        bb = []
        if hasfield(typeof(bright_bands), Symbol(Headers[k].direction))
            bb = getfield(bright_bands, Symbol(Headers[k].direction))
        end
        Headers[k].BeamData = synthesize_iq_data_vectrzd(;
            vr=vr,
            Header=H,
            spec_w=sv,
            SNR=SNR,
            add_iq_noise=add_iq_noise,
            dc_rel_power=dc_rel_power,
            gc_rel_power=gc_rel_power,
            sigma_dc=sigma_dc,
            bright_bands=bb
        )
    end

    if filepath != ""
        write_raw_file(Headers, filepath)
    end
    obj = Observation(Headers)

    obj.ref_height = z
    obj.ref_U = u
    obj.ref_V = v
    obj.ref_W = w

    return obj
end

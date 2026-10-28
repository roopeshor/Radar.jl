using .Radar: RadarData
using .io
# Observation  - value class representing all data in a single 5-beam radar observation.
#
# Each beam header is in north, east, west, south, vertical properties (all RadarData).
# Reference data "ref_..." comes from the .uvw file.
# Heights are from the .uvw file; moments references are per-beam.
# To use: instantiate with the base name of the observation files —
# the .raw, .uvw, and .mmts files are auto-located and parsed.
mutable struct Observation
    north::RadarData
    east::RadarData
    west::RadarData
    south::RadarData
    vertical::RadarData

    raw_filepath::String

    ref_height::Vector{Float64}
    ref_U::Vector{Float64}
    ref_V::Vector{Float64}
    ref_W::Vector{Float64}

    # computed wind components
    U::Vector{Float64}
    V::Vector{Float64}
    W::Vector{Float64}

    heights::Vector{Float64}     # from utils.compute_height_bins
    velocities::Vector{Float64}  # from utils.compute_velocity_bins

    # DBS factors [m/s per Hz]: λ/(2sin(oz)) and λ/2
    DBS_Factor_H::Float64
    DBS_Factor::Float64

    f_max::Float64        # maximum unambiguous Doppler frequency (Hz)
    v_max::Float64        # maximum unambiguous velocity (m/s)
    start_height::Float64
    end_height::Float64
    nRangeBins::Float64
    ipp_us::Float64
    in_coh::Float64
    n_coh::Float64
    nFFT::Float64

    function Observation()
        new(
            RadarData(), RadarData(), RadarData(), RadarData(), RadarData(),
            "",
            Float64[], Float64[], Float64[], Float64[],
            Float64[], Float64[], Float64[],
            Float64[], Float64[],
            0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
    end
end

# Outer constructor — mirrors MATLAB behaviour.
# `stuff` can be:
#   - a filepath String => reads .raw + attaches .mmts/.uvw references
#   - a Vector{RadarData} of length 5 => use directly (from synth, etc.)
function Observation(stuff=""; add_mmts::Bool=true, add_uvw::Bool=true)
    obs = Observation()

    if stuff == ""
        return obs
    end

    if stuff isa AbstractString
        obs.raw_filepath = endswith(stuff, ".raw") ? String(stuff) : String(stuff) * ".raw"

        if !isfile(obs.raw_filepath)
            error("File not found: $(obs.raw_filepath)")
        end

        raw      = io.read_raw_file(obs.raw_filepath)
        ref_data = io.add_reference_data(raw, obs.raw_filepath;
                       add_mmts=add_mmts, add_uvw=add_uvw)

        obs.north    = ref_data.beam_struct.north
        obs.east     = ref_data.beam_struct.east
        obs.west     = ref_data.beam_struct.west
        obs.south    = ref_data.beam_struct.south
        obs.vertical = ref_data.beam_struct.vertical

        obs.ref_height = ref_data.height
        obs.ref_U      = ref_data.U
        obs.ref_V      = ref_data.V
        obs.ref_W      = ref_data.W

    elseif stuff isa AbstractVector && eltype(stuff) <: RadarData && length(stuff) == 5
        for k in 1:5
            dir = stuff[k].direction
            if     dir == "north";    obs.north    = stuff[k]
            elseif dir == "east";     obs.east     = stuff[k]
            elseif dir == "west";     obs.west     = stuff[k]
            elseif dir == "south";    obs.south    = stuff[k]
            elseif dir == "vertical"; obs.vertical = stuff[k]
            end
        end
    else
        @warn "Observation: no suitable input data found"
        return obs
    end

    vel_factor       = 3.0e8 / 205e6 / 2.0
    obs.DBS_Factor   = vel_factor
    obs.DBS_Factor_H = vel_factor / sind(obs.north.m_fOffZenith)

    obs.f_max        = utils.compute_max_freq(obs.north.ipp_us, obs.north.n_coh)
    obs.v_max        = obs.f_max * vel_factor
    obs.start_height = obs.north.start_height
    obs.end_height   = obs.north.end_height
    obs.nRangeBins   = obs.north.nRangeBins
    obs.nFFT         = obs.north.nFFT
    obs.ipp_us       = obs.north.ipp_us
    obs.in_coh       = obs.north.in_coh
    obs.n_coh        = obs.north.n_coh

    obs.heights    = utils.compute_height_bins(obs.start_height, obs.end_height, obs.nRangeBins)
    obs.velocities = utils.compute_velocity_bins(obs.nFFT, obs.ipp_us, obs.n_coh)

    for dir in Data.directions
        beam_obj = getproperty(obs, Symbol(dir))
        beam_obj.spectra = utils.get_original_spectrum(beam_obj)
    end

    return obs
end

# Applies a moments/spectra function across all 5 beams in an Observation.
# Extra keyword args are forwarded to fx (mirrors MATLAB varargin).
function apply_beam_fx(obs::Observation, fx, args...; kwargs...)
    for dir in Data.directions
        sym = Symbol(dir)
        setfield!(obs, sym, fx(getproperty(obs, sym), args...; kwargs...))
    end
    return obs
end

# Computes U, V, W wind components from beam M1 moments.
function compute_uvw(obs::Observation)
    sum_M1 = obs.east.M1 .+ obs.west.M1 .+ obs.north.M1 .+ obs.south.M1
    cth    = cosd(obs.north.m_fOffZenith)

    obs.U  = .-(obs.east.M1  .- obs.west.M1)  .* obs.DBS_Factor_H
    obs.V  = .-(obs.north.M1 .- obs.south.M1) .* obs.DBS_Factor_H
    obs.W  = .-obs.DBS_Factor .* (cth .* sum_M1 .+ obs.vertical.M1) ./ (4 * cth^2 + 1)

    return obs
end

# RadarData - represents one beam's header and IQ data, plus post-processing fields.
#
# Uses COMPOSITION instead of inheritance (Julia has no concrete struct inheritance).
# All BeamHeader fields (m_fOffZenith, m_fAzimuth, etc.) are transparently forwarded
# via getproperty/setproperty!, so call sites are identical to MATLAB:  rd.m_fOffZenith

const _BEAM_HEADER_FIELDS = Set(fieldnames(BeamHeader))

mutable struct RadarData
    header::BeamHeader

    # ------------- Raw IQ data -------------------------------------------
    BeamData::Array{Complex{Float64}, 3}   # (RangeBins × NFFT × InCohIntegrations)

    # ---- aliases (mirrors of header fields, kept in sync by alias_params!) ----
    ipp_us::Float64        # m_fIntrPulsePeriod_us
    n_coh::Float64         # m_sNumOfCohIntegrations
    in_coh::Float64        # m_sNumOfInCohIntegrations
    start_height::Float64  # m_fWindow1StartHeight
    end_height::Float64    # m_fWindow1EndHeight
    nRangeBins::Float64    # m_sNumOfRangeBins
    nFFT::Float64          # m_sNFFT

    # --- Post-processing fields -------------------------------------------
    direction::String
    spectra::Matrix{Float64}          # (RangeBins × NFFT) averaged normalised power spectrum
    denoised_spectra::Matrix{Float64}
    M0::Vector{Float64}
    M1::Vector{Float64}
    M2::Vector{Float64}
    SNR::Vector{Float64}
    noise_level::Vector{Float64}

    a_params::DynamicProps   # algorithm-specific parameters — full dot-access notation

    # ---- reference data --------------------------------------------------
    ref_height::Vector{Float64}
    ref_M0::Vector{Float64}
    ref_M1::Vector{Float64}
    ref_M2::Vector{Float64}
    ref_SNR::Vector{Float64}
    ref_noise_level::Vector{Float64}

    function RadarData()
        rd = new()
        rd.header             = BeamHeader()
        rd.BeamData           = Array{Complex{Float64}}(undef, 0, 0, 0)
        rd.ipp_us             = 0.0
        rd.n_coh              = 0.0
        rd.in_coh             = 0.0
        rd.start_height       = 0.0
        rd.end_height         = 0.0
        rd.nRangeBins         = 0.0
        rd.nFFT               = 0.0
        rd.direction          = ""
        rd.spectra            = Matrix{Float64}(undef, 0, 0)
        rd.denoised_spectra   = Matrix{Float64}(undef, 0, 0)
        rd.M0                 = Float64[]
        rd.M1                 = Float64[]
        rd.M2                 = Float64[]
        rd.SNR                = Float64[]
        rd.noise_level        = Float64[]
        rd.a_params           = DynamicProps()
        rd.ref_height         = Float64[]
        rd.ref_M0             = Float64[]
        rd.ref_M1             = Float64[]
        rd.ref_M2             = Float64[]
        rd.ref_SNR            = Float64[]
        rd.ref_noise_level    = Float64[]
        return rd
    end
end

# --- Transparent property forwarding to BeamHeader ---
# Gives rd.m_fOffZenith syntax identical to MATLAB — no call-site changes needed.
function Base.getproperty(rd::RadarData, s::Symbol)
    s in _BEAM_HEADER_FIELDS && return getproperty(getfield(rd, :header), s)
    return getfield(rd, s)
end

function Base.setproperty!(rd::RadarData, s::Symbol, v)
    s in _BEAM_HEADER_FIELDS && return setproperty!(getfield(rd, :header), s, v)
    setfield!(rd, s, v)
end

# --- Constructor from any cfg (struct / NamedTuple / Dict) matching MATLAB behavior ---
function RadarData(cfg)
    rd = RadarData()

    # Populate BeamHeader fields
    for f in fieldnames(BeamHeader)
        if (cfg isa AbstractDict && haskey(cfg, String(f))) ||
           (cfg isa NamedTuple  && haskey(cfg, f)) ||
           (hasproperty(cfg, f))
            setproperty!(rd, f, getproperty(cfg, f))
        end
    end

    # Populate own fields (skip header and BeamData — handled separately)
    for f in fieldnames(RadarData)
        f in (:header, :BeamData, :a_params) && continue
        if (cfg isa AbstractDict && haskey(cfg, String(f))) ||
           (cfg isa NamedTuple  && haskey(cfg, f)) ||
           (hasproperty(cfg, f))
            setfield!(rd, f, getproperty(cfg, f))
        end
    end

    if (cfg isa AbstractDict && haskey(cfg, "BeamData")) ||
       (cfg isa NamedTuple  && haskey(cfg, :BeamData)) ||
       (hasproperty(cfg, :BeamData))
        rd.BeamData = getproperty(cfg, :BeamData)
    end

    alias_params!(rd)
    return rd
end

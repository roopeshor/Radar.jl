function create_synthetic_observation(
  u::AbstractVector,
  v::AbstractVector,
  w::AbstractVector;
  headerFields::AbstractVector{<:RadarData},
  filepath::AbstractString = "",
  sigma_v = nothing,
  SNR::AbstractVector = Float64[],
  sigma_dc::AbstractVector = Float64[],
  add_iq_noise::Bool = true,
  dc_rel_power::Real = 0.0,
  gc_rel_power::Real = 0.005,
  bright_bands = Dict{String, Matrix{Float64}}(),
  sigma_turb_base::Real = 0.0,
  sigma_chi::Real = 0.1,
)
  length(headerFields) == 5 ||
    throw(ArgumentError("headerFields must contain five beams"))
  headers = deepcopy(collect(headerFields))
  first_header = headers[1]
  range_res = 3e8 * first_header.m_fBaudLength_us * 1e-6 / 2
  z =
    first_header.m_fWindow1StartHeight .+
    collect(0:(Int(first_header.m_sNumOfRangeBins) - 1)) .* range_res
  dbs_v = 3e8/205e6/2

  for k in eachindex(headers)
    h = headers[k]
    vr = Utils.compute_radial_velocity(u, v, w, h.m_fAzimuth, h.m_fOffZenith)
    sw = sigma_v
    if sw === nothing || sw isa AbstractString
      sw = compute_spectrum_broadening(
        z,
        u,
        v,
        w,
        vr;
        sigma_turb_base = sigma_turb_base,
        sigma_chi = sigma_chi,
      )
    end
    h.ref_M1 = -Float64.(vr) ./ dbs_v
    bb = get(bright_bands, h.direction, zeros(0, 2))
    h.BeamData = synthesize_iq_data(
      vr = vr,
      Header = h,
      spec_w = sw,
      SNR = SNR,
      add_iq_noise = add_iq_noise,
      dc_rel_power = dc_rel_power,
      gc_rel_power = gc_rel_power,
      sigma_dc = sigma_dc,
      bright_bands = bb,
    )
  end
  isempty(filepath) || write_raw_file(headers, filepath)

  obs = Observation()
  for h in headers
    setproperty!(obs, Symbol(h.direction), h)
  end
  obs.ref_height = z
  obs.ref_U = Float64.(u)
  obs.ref_V = Float64.(v)
  obs.ref_W = Float64.(w)
  obs.raw_filepath = String(filepath)
  vh = obs.north
  factor = 3e8/205e6/2
  obs.DBS_Factor = factor
  obs.DBS_Factor_H = factor/(2sind(vh.m_fOffZenith))
  obs.f_max = Utils.compute_max_freq(vh.ipp_us, vh.n_coh)
  obs.v_max = obs.f_max*factor
  obs.ipp_us = vh.m_fIntrPulsePeriod_us
  obs.n_coh = vh.m_sNumOfCohIntegrations
  obs.in_coh = vh.m_sNumOfInCohIntegrations
  obs.start_height = vh.m_fWindow1StartHeight
  obs.end_height = vh.m_fWindow1EndHeight
  obs.nRangeBins = vh.m_sNumOfRangeBins
  obs.nFFT = vh.m_sNFFT
  obs.heights = Utils.compute_height_bins(
    obs.start_height,
    obs.end_height,
    Int(obs.nRangeBins),
  )
  obs.velocities =
    Utils.compute_velocity_bins(Int(obs.nFFT), obs.ipp_us, obs.n_coh)
  for direction in Data.directions
    beam = getproperty(obs, direction)
    beam.spectra = Utils.get_original_spectrum(beam)
  end
  return obs
end

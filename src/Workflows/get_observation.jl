function _populate_observation!(
  obs::Observation,
  beams::AbstractVector{<:RadarData};
  raw_filepath::AbstractString = "",
)::Observation
  if length(beams) != 5
    throw(ArgumentError("Expected five radar beams, got $(length(beams))"))
  end
  for beam in beams
    direction = beam.direction
    if !(Symbol(direction) in Data.directions)
      throw(ArgumentError("Unrecognized beam direction '$direction'"))
    end
    setproperty!(obs, Symbol(direction), beam)
  end
  obs.raw_filepath = String(raw_filepath)
  north = obs.north
  factor = 3.0e8/205e6/2

  obs.DBS_Factor   = factor
  obs.DBS_Factor_H = factor/(2*sind(north.m_fOffZenith))
  obs.f_max        = Utils.compute_max_freq(north.ipp_us, north.n_coh)
  obs.v_max        = obs.f_max*factor
  obs.ipp_us       = north.m_fIntrPulsePeriod_us
  obs.n_coh        = north.m_sNumOfCohIntegrations
  obs.in_coh       = north.m_sNumOfInCohIntegrations
  obs.start_height = north.m_fWindow1StartHeight
  obs.end_height   = north.m_fWindow1EndHeight
  obs.nRangeBins   = north.m_sNumOfRangeBins
  obs.nFFT         = north.m_sNFFT

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

get_observation() = Observation()

function get_observation(
  stuff::AbstractString;
  add_mmts::Bool = true,
  add_uvw::Bool = true,
)::Observation
  isempty(stuff) && return Observation()
  raw_path = endswith(stuff, ".raw") ? String(stuff) : String(stuff)*".raw"
  beams = RadarIO.read_raw_file(raw_path)
  isempty(beams) && throw(ArgumentError("No beams were read from $raw_path"))
  refs = RadarIO.add_reference_data(beams, raw_path; add_mmts, add_uvw)
  obs = Observation()
  referenced_beams = RadarData[refs.beam_struct[d] for d in Data.directions]

  obs.ref_height = refs.height
  obs.ref_U      = refs.U
  obs.ref_V      = refs.V
  obs.ref_W      = refs.W
  return _populate_observation!(obs, referenced_beams; raw_filepath = raw_path)
end

function get_observation(beams::AbstractVector{<:RadarData})::Observation
  obs = Observation()
  return _populate_observation!(obs, beams)
end

# Preserve the original MATLAB-facing spelling as a compatibility alias.
getObservation(args...; kwargs...) = get_observation(args...; kwargs...)

"Attach `.mmts` beam moments and `.uvw` wind reference data to a beam collection."
function add_reference_data(
  beams::AbstractVector{<:RadarData},
  filepath::AbstractString;
  add_mmts::Bool = true,
  add_uvw::Bool = true,
  uvw_file_ext::AbstractString = "_W1.uvw",
  mmts_file_ext::AbstractDict = Data.dir2mmts_file_ext,
)
  length(beams) == 5 ||
    throw(ArgumentError("Expected five beams, received $(length(beams))"))
  folder = dirname(filepath)
  base = splitext(basename(filepath))[1]
  beam_struct = Dict{Symbol, RadarData}(
    :vertical => RadarData(),
    :north    => RadarData(),
    :west     => RadarData(),
    :south    => RadarData(),
    :east     => RadarData(),
  )

  ref_height = Float64[]
  U = Float64[]
  V = Float64[]
  W = Float64[]
  uvw_path = joinpath(folder, base * uvw_file_ext)
  for source_beam in beams
    beam = deepcopy(source_beam)
    direction = beam.direction

    beam.ref_height      = Float64[]
    beam.ref_M0          = Float64[]
    beam.ref_M1          = Float64[]
    beam.ref_M2          = Float64[]
    beam.ref_SNR         = Float64[]
    beam.ref_noise_level = Float64[]

    if add_mmts
      suffix = get(mmts_file_ext, direction, nothing)
      suffix === nothing && throw(
        ArgumentError("No .mmts suffix configured for direction '$direction'"),
      )
      path = joinpath(folder, base * suffix)
      if isfile(path)
        headers, table = _read_numeric_table(path)

        beam.direction       = Data.aoz2dir(az = beam.m_fAzimuth, oz = beam.m_fOffZenith)
        beam.ref_height      = _column(headers, table, "Height (km)")
        beam.ref_M0          = _column(headers, table, "M0 Total Power(dBm)")
        beam.ref_M1          = _column(headers, table, "M1 Mean Doppler (Hz)")
        beam.ref_M2          = _column(headers, table, "M2 Doppler Spread (Hz)")
        beam.ref_SNR         = _column(headers, table, "SNR (dB)")
        beam.ref_noise_level = _column(headers, table, "Noise Level (dBm)")
      else
        @warn "Reference moment file not found; reference arrays left empty" path
      end
    end
    beam_struct[Symbol(direction)] = beam
  end

  if add_uvw
    if isfile(uvw_path)
      headers, table = _read_numeric_table(uvw_path)
      ref_height = _column(headers, table, "Height (km)")
      U = _column(headers, table, "Zonal(U) (m/s)")
      V = _column(headers, table, "Meridional(V) (m/s)")
      W = _column(headers, table, "Vertical (W) (m/s)")
    else
      @warn "UVW reference file not found; reference wind arrays left empty" uvw_path
    end
  end
  return (; beam_struct, height = ref_height, U, V, W)
end

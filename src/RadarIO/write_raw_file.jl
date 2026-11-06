"Write a vector of beams to the instrument's 1024-byte-header `.raw` format."
function write_raw_file(
  beams::AbstractVector{<:RadarData},
  filepath::AbstractString,
)
  isempty(beams) && throw(ArgumentError("No radar beams supplied"))
  dir = dirname(filepath)
  isempty(dir) || mkpath(dir)
  open(filepath, "w") do io
    nbeams = length(beams)
    for (beam_no, b) in enumerate(beams)
      ni = Int(b.m_sNumOfInCohIntegrations)
      ni <= 0 && (ni = 1)
      nr = Int(b.m_sNumOfRangeBins)
      nf = Int(b.m_sNFFT)
      if !isempty(b.BeamData)
        size(b.BeamData, 1) == nr && size(b.BeamData, 2) == nf || throw(
          DimensionMismatch(
            "BeamData dimensions do not match header dimensions",
          ),
        )
        size(b.BeamData, 3) >= ni || throw(
          DimensionMismatch(
            "BeamData has fewer incoherent integrations than the header",
          ),
        )
      end
      for integration = 1:ni
        h = zeros(UInt8, 1024)
        int_fields = [
          (1, b.m_sMagicNumber),
          (3, b.m_sNumOfRangeBins),
          (9, b.m_sNFFT),
          (11, b.m_sNumOfCohIntegrations),
          (13, b.m_sNumOfInCohIntegrations),
          (15, b.m_sCodeFlag),
          (25, b.m_sNumOfObservWindows),
          (27, b.m_sStationID),
          (29, b.m_sReserved1),
          (31, b.m_sReserved2),
          (33, b.m_sReserved3),
          (35, b.m_sYear),
          (37, b.m_sMonth),
          (39, b.m_sDay),
          (41, b.m_sHour),
          (43, b.m_sMin),
          (45, b.m_sSec),
          (47, b.m_sReserved4),
          (117, b.m_sOperationMode),
          # The source MATLAB writer writes the 1-based beam ordinal here.
          (119, beam_no),
          (121, nbeams),
          (123, b.m_usCodeLength),
          (149, b.m_usWindowType),
          (151, b.m_usReserved),
        ]
        for (offset, value) in int_fields
          _putle!(h, offset, Int16(round(value)))
        end
        float_fields = [
          (5, b.m_fBaudLength_us),
          (17, b.m_fIntrPulsePeriod_us),
          (21, b.m_fPulseWidth_us),
          (49, b.m_fLatitude),
          (53, b.m_fLongitude),
          (57, b.m_fAltitude),
          (61, b.m_fTotalPowerRadiated),
          (141, b.m_fTRP_RF_Delay_us),
          (145, b.m_fMGC_Gain_dB),
          (685, b.m_fAzimuth),
          (689, b.m_fOffZenith),
          (693, b.m_fWindow1StartHeight),
          (697, b.m_fWindow1EndHeight),
          (701, b.m_fWindow2StartHeight),
          (705, b.m_fWindow2EndHeight),
          (709, b.m_fWindow3StartHeight),
          (713, b.m_fWindow3EndHeight),
          (717, b.m_fWindow4StartHeight),
          (721, b.m_fWindow4EndHeight),
          (725, b.m_fWindow5StartHeight),
          (729, b.m_fWindow5EndHeight),
        ]
        for (offset, value) in float_fields
          _putle!(h, offset, Float32(value))
        end
        pwr =
          isempty(b.m_fRadiatedPower_ClusterWise) ? zeros(13) :
          b.m_fRadiatedPower_ClusterWise
        _put_vector!(h, 65, pwr, Float32)
        copyto!(h, 125, _fixed_bytes(b.m_cComment, 16), 1, 16)
        copyto!(h, 153, _fixed_bytes(b.m_cArrRemarks, 512), 1, 512)
        _put_vector!(
          h,
          665,
          isempty(b.m_ulCodeA) ? zeros(UInt32, 2) : b.m_ulCodeA,
          UInt32,
        )
        _put_vector!(
          h,
          673,
          isempty(b.m_ulCodeB) ? zeros(UInt32, 2) : b.m_ulCodeB,
          UInt32,
        )
        _putle!(h, 681, Int8(round(b.m_ucExperimentDataEnable)))
        _put_vector!(
          h,
          682,
          isempty(b.m_cReserved) ? zeros(Int8, 3) : b.m_cReserved,
          Int8,
        )
        write(io, h)

        if isempty(b.BeamData)
          rawzero = zeros(UInt8, 8 * nr * nf)
          write(io, rawzero)
        else
          block = vec(@view b.BeamData[:, :, integration])
          raw = Vector{Int32}(undef, 2 * length(block))
          for i in eachindex(block)
            raw[2i - 1] = _round_away_int32(real(block[i]))
            raw[2i] = _round_away_int32(imag(block[i]))
          end
          rawbytes = zeros(UInt8, 8 * length(block))
          for i in eachindex(raw)
            _putle!(rawbytes, 4i-3, raw[i])
          end
          write(io, rawbytes)
        end
      end
    end
  end
  return nothing
end

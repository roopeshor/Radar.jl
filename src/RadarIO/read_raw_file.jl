"Read a radar `.raw` file into a vector of `RadarData` beams."
function read_raw_file(filepath::AbstractString)
  isfile(filepath) || throw(ArgumentError("Cannot open file: $filepath"))
  beams = RadarData[]
  open(filepath, "r") do io
    seek(io, 120) # zero-based byte offset 120: beam count at bytes 121–122
    count_bytes = read(io, 2)
    length(count_bytes) == 2 || return beams
    beam_count = Int(_readle(Int16, count_bytes, 1))
    beam_count <= 0 && return beams
    seekstart(io)

    for beam_index = 1:beam_count
      header = read(io, 1024)
      length(header) == 1024 || break
      magic = Int(_readle(Int16, header, 1))
      magic == 369 || throw(
        ArgumentError(
          "Unexpected radar header magic number $magic in beam $beam_index",
        ),
      )

      b = RadarData()
      b.m_sMagicNumber = magic
      # Header fields: each offset is one-based, matching the MATLAB layout.
      int16_fields = [
        (:m_sNumOfRangeBins, 3),
        (:m_sNFFT, 9),
        (:m_sNumOfCohIntegrations, 11),
        (:m_sNumOfInCohIntegrations, 13),
        (:m_sCodeFlag, 15),
        (:m_sNumOfObservWindows, 25),
        (:m_sStationID, 27),
        (:m_sReserved1, 29),
        (:m_sReserved2, 31),
        (:m_sReserved3, 33),
        (:m_sYear, 35),
        (:m_sMonth, 37),
        (:m_sDay, 39),
        (:m_sHour, 41),
        (:m_sMin, 43),
        (:m_sSec, 45),
        (:m_sReserved4, 47),
        (:m_sOperationMode, 117),
        (:m_sCurrentBeamCnt, 119),
        (:m_sTotalNumberofBeams, 121),
        (:m_usCodeLength, 123),
        (:m_usWindowType, 149),
        (:m_usReserved, 151),
      ]
      for (name, offset) in int16_fields
        setproperty!(b, name, Float64(_readle(Int16, header, offset)))
      end
      float32_fields = [
        (:m_fBaudLength_us, 5),
        (:m_fIntrPulsePeriod_us, 17),
        (:m_fPulseWidth_us, 21),
        (:m_fLatitude, 49),
        (:m_fLongitude, 53),
        (:m_fAltitude, 57),
        (:m_fTotalPowerRadiated, 61),
        (:m_fTRP_RF_Delay_us, 141),
        (:m_fMGC_Gain_dB, 145),
        (:m_fAzimuth, 685),
        (:m_fOffZenith, 689),
        (:m_fWindow1StartHeight, 693),
        (:m_fWindow1EndHeight, 697),
        (:m_fWindow2StartHeight, 701),
        (:m_fWindow2EndHeight, 705),
        (:m_fWindow3StartHeight, 709),
        (:m_fWindow3EndHeight, 713),
        (:m_fWindow4StartHeight, 717),
        (:m_fWindow4EndHeight, 721),
        (:m_fWindow5StartHeight, 725),
        (:m_fWindow5EndHeight, 729),
      ]
      for (name, offset) in float32_fields
        setproperty!(b, name, Float64(_readle(Float32, header, offset)))
      end
      b.m_fRadiatedPower_ClusterWise =
        Float64[_readle(Float32, header, i) for i = 65:4:113]
      b.m_cComment = Vector{UInt8}(header[125:140])
      b.m_cArrRemarks = Vector{UInt8}(header[153:664])
      b.m_ulCodeA =
        UInt32[_readle(UInt32, header, 665), _readle(UInt32, header, 669)]
      b.m_ulCodeB =
        UInt32[_readle(UInt32, header, 673), _readle(UInt32, header, 677)]
      b.m_ucExperimentDataEnable = Float64(_readle(Int8, header, 681))
      b.m_cReserved = Int8[
        _readle(Int8, header, 682),
        _readle(Int8, header, 683),
        _readle(Int8, header, 684),
      ]

      b.m_sNumOfInCohIntegrations <= 0 && (b.m_sNumOfInCohIntegrations = 1)
      nr = Int(b.m_sNumOfRangeBins)
      nf = Int(b.m_sNFFT)
      ni = Int(b.m_sNumOfInCohIntegrations)
      nr > 0 && nf > 0 && ni > 0 || throw(
        ArgumentError(
          "Invalid dimensions in beam $beam_index: ($nr, $nf, $ni)",
        ),
      )
      b.BeamData = Array{ComplexF64}(undef, nr, nf, ni)
      add_extra_RD_params!(b)

      for integration = 1:ni
        if integration > 1
          extra_header = read(io, 1024)
          length(extra_header) == 1024 || throw(EOFError())
        end
        rawbytes = read(io, 8 * nr * nf)
        length(rawbytes) == 8 * nr * nf || throw(EOFError())
        raw = Int32[_readle(Int32, rawbytes, 4i-3) for i = 1:(2 * nr * nf)]
        iq = ComplexF64.(raw[1:2:end], raw[2:2:end])
        b.BeamData[:, :, integration] .= reshape(iq, nr, nf)
      end
      push!(beams, b)
    end
  end
  return beams
end

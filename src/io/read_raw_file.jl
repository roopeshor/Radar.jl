using ..Radar: RadarData

# io.read_raw_file Binary parser for ST radar (.raw) experimental data files.
#
#   Reads 1024-byte binary header blocks and 32-bit complex IQ time series data for all beams,
#   constructing an array of RadarData objects populated with beam metadata and BeamData cubes.
function read_raw_file(filepath::String)::NTuple{5, RadarData}
    if !isfile(filepath)
        error("Cannot open file: $filepath")
    end

    fPtr = open(filepath, "r")
    DEFAULT_MAGIC_NUMBER = 369

    seek(fPtr, 120)
    beam_count = Int(read(fPtr, Int16))
    seek(fPtr, 0)

    if beam_count <= 0
        close(fPtr)
        return RadarData[]
    end

    Headers = RadarData[]

    for beam_No in 1:beam_count
        hdrBytes = read(fPtr, 1024)
        if length(hdrBytes) < 1024
            break
        end

        magicNum = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 1))))

        rd = RadarData()
        rd.m_sMagicNumber = magicNum

        if magicNum != DEFAULT_MAGIC_NUMBER
            @error "Wrong File Selected: ST Radar Warning"
            close(fPtr)
            return RadarData[]
        end

        rd.m_sNumOfRangeBins            = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 3))))
        rd.m_fBaudLength_us             = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 5))))
        rd.m_sNFFT                      = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 9))))
        rd.m_sNumOfCohIntegrations      = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 11))))
        rd.m_sNumOfInCohIntegrations    = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 13))))
        rd.m_sCodeFlag                  = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 15))))
        rd.m_fIntrPulsePeriod_us        = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 17))))
        rd.m_fPulseWidth_us             = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 21))))
        rd.m_sNumOfObservWindows        = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 25))))
        rd.m_sStationID                 = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 27))))
        rd.m_sReserved1                 = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 29))))
        rd.m_sReserved2                 = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 31))))
        rd.m_sReserved3                 = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 33))))
        rd.m_sYear                      = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 35))))
        rd.m_sMonth                     = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 37))))
        rd.m_sDay                       = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 39))))
        rd.m_sHour                      = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 41))))
        rd.m_sMin                       = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 43))))
        rd.m_sSec                       = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 45))))
        rd.m_sReserved4                 = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 47))))
        rd.m_fLatitude                  = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 49))))
        rd.m_fLongitude                 = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 53))))
        rd.m_fAltitude                  = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 57))))
        rd.m_fTotalPowerRadiated        = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 61))))

        cluster_ptr = reinterpret(Ptr{Float32}, pointer(hdrBytes, 65))
        rd.m_fRadiatedPower_ClusterWise = [Float64(unsafe_load(cluster_ptr, i)) for i in 1:13]

        rd.m_sOperationMode             = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 117))))
        rd.m_sCurrentBeamCnt            = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 119))))
        rd.m_sTotalNumberofBeams        = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 121))))
        rd.m_usCodeLength               = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 123))))

        comment_ptr = reinterpret(Ptr{Int8}, pointer(hdrBytes, 125))
        rd.m_cComment                   = [Char(unsafe_load(comment_ptr, i)) for i in 1:16]

        rd.m_fTRP_RF_Delay_us           = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 141))))
        rd.m_fMGC_Gain_dB               = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 145))))
        rd.m_usWindowType               = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 149))))
        rd.m_usReserved                 = Float64(unsafe_load(reinterpret(Ptr{Int16}, pointer(hdrBytes, 151))))

        remarks_ptr = reinterpret(Ptr{Int8}, pointer(hdrBytes, 153))
        rd.m_cArrRemarks                = [Char(unsafe_load(remarks_ptr, i)) for i in 1:512]

        codea_ptr = reinterpret(Ptr{UInt32}, pointer(hdrBytes, 665))
        rd.m_ulCodeA                    = [Float64(unsafe_load(codea_ptr, i)) for i in 1:2]

        codeb_ptr = reinterpret(Ptr{UInt32}, pointer(hdrBytes, 673))
        rd.m_ulCodeB                    = [Float64(unsafe_load(codeb_ptr, i)) for i in 1:2]

        rd.m_ucExperimentDataEnable     = Float64(unsafe_load(reinterpret(Ptr{Int8}, pointer(hdrBytes, 681))))

        res_bytes_ptr = reinterpret(Ptr{Int8}, pointer(hdrBytes, 682))
        rd.m_cReserved                  = [unsafe_load(res_bytes_ptr, i) for i in 1:3]

        rd.m_fAzimuth                   = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 685))))
        rd.m_fOffZenith                 = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 689))))
        rd.m_fWindow1StartHeight        = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 693))))
        rd.m_fWindow1EndHeight          = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 697))))
        rd.m_fWindow2StartHeight        = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 701))))
        rd.m_fWindow2EndHeight          = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 705))))
        rd.m_fWindow3StartHeight        = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 709))))
        rd.m_fWindow3EndHeight          = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 713))))
        rd.m_fWindow4StartHeight        = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 717))))
        rd.m_fWindow4EndHeight          = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 721))))
        rd.m_fWindow5StartHeight        = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 725))))
        rd.m_fWindow5EndHeight          = Float64(unsafe_load(reinterpret(Ptr{Float32}, pointer(hdrBytes, 729))))

        if rd.m_sNumOfInCohIntegrations <= 0
            rd.m_sNumOfInCohIntegrations = 1.0
        end

        nRangeBins = Int(rd.m_sNumOfRangeBins)
        nNFFT = Int(rd.m_sNFFT)
        nInCoh = Int(rd.m_sNumOfInCohIntegrations)
        dataBlockSize = nRangeBins * nNFFT

        rd.BeamData = Array{Complex{Float64}}(undef, nRangeBins, nNFFT, nInCoh)

        for InCohIntegration in 1:nInCoh
            if InCohIntegration != 1
                headerBytesData = read(fPtr, 1024)
                if length(headerBytesData) != 1024
                    @error "Error in file Read: Expected Header Data Count 1024"
                end
            end

            raw_data_vector = reinterpret(Int32, read(fPtr, 2 * dataBlockSize * sizeof(Int32)))
            rawData = reshape(raw_data_vector, 2, dataBlockSize)

            in_data = complex.(Float64.(rawData[1, :]), Float64.(rawData[2, :]))
            rd.BeamData[:, :, InCohIntegration] = reshape(in_data, (nRangeBins, nNFFT))
        end
		rdd = RadarData(rd);
		push!(Headers, rdd);
    end

    close(fPtr)
    return Tuple(Headers)::NTuple{5, RadarData}
end

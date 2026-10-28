# Binary serializer for ST radar (.raw) experimental data files.
#
#   Serializes an array of RadarData objects into standard 1024-byte binary headers and interleaved
#   32-bit real/imaginary complex IQ time series data blocks.
function write_raw_file(Header::Vector{RadarData}, filepath::String)
    fPtr = open(filepath, "w")

    beam_count = length(Header)

    for beam_No in 1:beam_count
        nInCoh = Int(Header[beam_No].m_sNumOfInCohIntegrations)
        if nInCoh <= 0
            nInCoh = 1
        end

        for InCohIntegration in 1:nInCoh
            hdrBytes = zeros(UInt8, 1024)

            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sMagicNumber), 1, 2)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sNumOfRangeBins), 3, 4)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fBaudLength_us), 5, 8)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sNFFT), 9, 10)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sNumOfCohIntegrations), 11, 12)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sNumOfInCohIntegrations), 13, 14)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sCodeFlag), 15, 16)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fIntrPulsePeriod_us), 17, 20)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fPulseWidth_us), 21, 24)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sNumOfObservWindows), 25, 26)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sStationID), 27, 28)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sReserved1), 29, 30)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sReserved2), 31, 32)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sReserved3), 33, 34)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sYear), 35, 36)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sMonth), 37, 38)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sDay), 39, 40)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sHour), 41, 42)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sMin), 43, 44)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sSec), 45, 46)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sReserved4), 47, 48)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fLatitude), 49, 52)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fLongitude), 53, 56)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fAltitude), 57, 60)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fTotalPowerRadiated), 61, 64)

            pow_cluster = Header[beam_No].m_fRadiatedPower_ClusterWise
            if isempty(pow_cluster)
                pow_cluster = zeros(Float32, 13)
            else
                pow_cluster = Float32.(pow_cluster)
            end
            hdrBytes = assign_bytes(hdrBytes, pow_cluster, 65, 116)

            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_sOperationMode), 117, 118)
            hdrBytes = assign_bytes(hdrBytes, Int16(beam_No), 119, 120)
            hdrBytes = assign_bytes(hdrBytes, Int16(beam_count), 121, 122)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_usCodeLength), 123, 124)

            cmt = Header[beam_No].m_cComment
            if isempty(cmt)
                cmt_vec = zeros(Int8, 16)
            else
                cmt_vec = Int8.([c for c in cmt])
                if length(cmt_vec) < 16
                    append!(cmt_vec, zeros(Int8, 16 - length(cmt_vec)))
                else
                    cmt_vec = cmt_vec[1:16]
                end
            end
            hdrBytes = assign_bytes(hdrBytes, cmt_vec, 125, 140)

            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fTRP_RF_Delay_us), 141, 144)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fMGC_Gain_dB), 145, 148)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_usWindowType), 149, 150)
            hdrBytes = assign_bytes(hdrBytes, Int16(Header[beam_No].m_usReserved), 151, 152)

            rem = Header[beam_No].m_cArrRemarks
            if isempty(rem)
                rem_vec = zeros(Int8, 512)
            else
                rem_vec = Int8.([c for c in rem])
                if length(rem_vec) < 512
                    append!(rem_vec, zeros(Int8, 512 - length(rem_vec)))
                else
                    rem_vec = rem_vec[1:512]
                end
            end
            hdrBytes = assign_bytes(hdrBytes, rem_vec, 153, 664)

            coda = Header[beam_No].m_ulCodeA
            if isempty(coda)
                coda = zeros(UInt32, 2)
            else
                coda = UInt32.(coda)
            end
            hdrBytes = assign_bytes(hdrBytes, coda, 665, 672)

            codb = Header[beam_No].m_ulCodeB
            if isempty(codb)
                codb = zeros(UInt32, 2)
            else
                codb = UInt32.(codb)
            end
            hdrBytes = assign_bytes(hdrBytes, codb, 673, 680)

            hdrBytes = assign_bytes(hdrBytes, Int8(Header[beam_No].m_ucExperimentDataEnable), 681, 681)

            cres = Header[beam_No].m_cReserved
            if isempty(cres)
                cres = zeros(Int8, 3)
            else
                cres = Int8.(cres)
            end
            hdrBytes = assign_bytes(hdrBytes, cres, 682, 684)

            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fAzimuth), 685, 688)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fOffZenith), 689, 692)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow1StartHeight), 693, 696)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow1EndHeight), 697, 700)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow2StartHeight), 701, 704)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow2EndHeight), 705, 708)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow3StartHeight), 709, 712)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow3EndHeight), 713, 716)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow4StartHeight), 717, 720)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow4EndHeight), 721, 724)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow5StartHeight), 725, 728)
            hdrBytes = assign_bytes(hdrBytes, Float32(Header[beam_No].m_fWindow5EndHeight), 729, 732)

            write(fPtr, hdrBytes)

            nRangeBins = Int(Header[beam_No].m_sNumOfRangeBins)
            nNFFT = Int(Header[beam_No].m_sNFFT)

            if !isempty(Header[beam_No].BeamData)
                block = Header[beam_No].BeamData[:, :, InCohIntegration]
                flat_block = vec(block)
                out_data = Array{Int32}(undef, 2 * length(flat_block))
                out_data[1:2:end] .= Int32.(real.(flat_block))
                out_data[2:2:end] .= Int32.(imag.(flat_block))

                write(fPtr, out_data)
            else
                write(fPtr, zeros(Int32, 2 * nRangeBins * nNFFT))
            end
        end
    end
    close(fPtr)
end

function assign_bytes(bytes::Vector{UInt8}, data, start_idx::Int, end_idx::Int)
    raw_bytes = reinterpret(UInt8, pointer_from_objref(Ref(data)))
    # For safety with primitive/array types:
    raw_bytes_vec = Vector{UInt8}(undef, end_idx - start_idx + 1)

    # Generic serialization into byte chunks
    io_buf = IOBuffer()
    write(io_buf, data)
    seekstart(io_buf)
    read!(io_buf, raw_bytes_vec)

    bytes[start_idx:end_idx] .= raw_bytes_vec
    return bytes
end

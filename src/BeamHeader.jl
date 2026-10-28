# BeamHeader - all 52 basic beam header data fields
# Plain concrete struct — no fields changed, just used via composition in RadarData
mutable struct BeamHeader
    m_sMagicNumber::Float64
    m_sNumOfRangeBins::Float64
    m_fBaudLength_us::Float64
    m_sNFFT::Float64
    m_sNumOfCohIntegrations::Float64
    m_sNumOfInCohIntegrations::Float64
    m_sCodeFlag::Float64
    m_fIntrPulsePeriod_us::Float64
    m_fPulseWidth_us::Float64
    m_sNumOfObservWindows::Float64
    m_sStationID::Float64
    m_sReserved1::Float64
    m_sReserved2::Float64
    m_sReserved3::Float64
    m_sYear::Float64
    m_sMonth::Float64
    m_sDay::Float64
    m_sHour::Float64
    m_sMin::Float64
    m_sSec::Float64
    m_sReserved4::Float64
    m_fLatitude::Float64
    m_fLongitude::Float64
    m_fAltitude::Float64
    m_fTotalPowerRadiated::Float64
    m_fRadiatedPower_ClusterWise::Vector{Float64}
    m_sOperationMode::Float64
    m_sCurrentBeamCnt::Float64
    m_sTotalNumberofBeams::Float64
    m_usCodeLength::Float64
    m_cComment::Vector{Char}
    m_fTRP_RF_Delay_us::Float64
    m_fMGC_Gain_dB::Float64
    m_usWindowType::Float64
    m_usReserved::Float64
    m_cArrRemarks::Vector{Char}
    m_ulCodeA::Vector{Float64}
    m_ulCodeB::Vector{Float64}
    m_ucExperimentDataEnable::Float64
    m_cReserved::Vector{Int8}
    m_fAzimuth::Float64
    m_fOffZenith::Float64
    m_fWindow1StartHeight::Float64
    m_fWindow1EndHeight::Float64
    m_fWindow2StartHeight::Float64
    m_fWindow2EndHeight::Float64
    m_fWindow3StartHeight::Float64
    m_fWindow3EndHeight::Float64
    m_fWindow4StartHeight::Float64
    m_fWindow4EndHeight::Float64
    m_fWindow5StartHeight::Float64
    m_fWindow5EndHeight::Float64

    function BeamHeader()
        new(
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, zeros(Float64, 13), 0.0, 0.0, 0.0, 0.0,
            fill(' ', 16), 0.0, 0.0, 0.0, 0.0, fill(' ', 512), zeros(Float64, 2), zeros(Float64, 2), 0.0,
            zeros(Int8, 3), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
    end
end

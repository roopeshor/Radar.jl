"""
Mutable representation of one ST-radar beam.

Header scalar fields are stored as Float64 to mirror the MATLAB class. Binary
serialization converts them back to their on-disk integer/Float32 widths.
Complex IQ samples are stored as ComplexF64 for consistency with MATLAB's
`double` complex arrays.
"""
Base.@kwdef mutable struct RadarData
  # Basic 1024-byte header fields
  m_sMagicNumber::Float64                       = 0
  m_sNumOfRangeBins::Float64                    = 0
  m_fBaudLength_us::Float64                     = 0
  m_sNFFT::Float64                              = 1024
  m_sNumOfCohIntegrations::Float64              = 0
  m_sNumOfInCohIntegrations::Float64            = 1
  m_sCodeFlag::Float64                          = 0
  m_fIntrPulsePeriod_us::Float64                = 0
  m_fPulseWidth_us::Float64                     = 0
  m_sNumOfObservWindows::Float64                = 0
  m_sStationID::Float64                         = 0
  m_sReserved1::Float64                         = 0
  m_sReserved2::Float64                         = 0
  m_sReserved3::Float64                         = 0
  m_sYear::Float64                              = 0
  m_sMonth::Float64                             = 0
  m_sDay::Float64                               = 0
  m_sHour::Float64                              = 0
  m_sMin::Float64                               = 0
  m_sSec::Float64                               = 0
  m_sReserved4::Float64                         = 0
  m_fLatitude::Float64                          = 0
  m_fLongitude::Float64                         = 0
  m_fAltitude::Float64                          = 0
  m_fTotalPowerRadiated::Float64                = 0
  m_fRadiatedPower_ClusterWise::Vector{Float64} = zeros(13)
  m_sOperationMode::Float64                     = 0
  m_sCurrentBeamCnt::Float64                    = 0
  m_sTotalNumberofBeams::Float64                = 0
  m_usCodeLength::Float64                       = 0
  m_cComment::Vector{UInt8}                     = zeros(UInt8, 16)
  m_fTRP_RF_Delay_us::Float64                   = 0
  m_fMGC_Gain_dB::Float64                       = 0
  m_usWindowType::Float64                       = 0
  m_usReserved::Float64                         = 0
  m_cArrRemarks::Vector{UInt8}                  = zeros(UInt8, 512)
  m_ulCodeA::Vector{UInt32}                     = zeros(UInt32, 2)
  m_ulCodeB::Vector{UInt32}                     = zeros(UInt32, 2)
  m_ucExperimentDataEnable::Float64             = 0
  m_cReserved::Vector{Int8}                     = zeros(Int8, 3)
  m_fAzimuth::Float64                           = 0
  m_fOffZenith::Float64                         = 0
  m_fWindow1StartHeight::Float64                = 0
  m_fWindow1EndHeight::Float64                  = 0
  m_fWindow2StartHeight::Float64                = 0
  m_fWindow2EndHeight::Float64                  = 0
  m_fWindow3StartHeight::Float64                = 0
  m_fWindow3EndHeight::Float64                  = 0
  m_fWindow4StartHeight::Float64                = 0
  m_fWindow4EndHeight::Float64                  = 0
  m_fWindow5StartHeight::Float64                = 0
  m_fWindow5EndHeight::Float64                  = 0

  # Raw IQ cube and convenience aliases
  BeamData::Array{ComplexF64, 3} = Array{ComplexF64}(undef, 0, 0, 0)
  ipp_us::Float64                = 0
  n_coh::Float64                 = 0
  in_coh::Float64                = 1
  start_height::Float64          = 0
  end_height::Float64            = 0
  nRangeBins::Float64            = 0
  nFFT::Float64                  = 1024

  # Processing fields
  direction::String                 = "unknown"
  spectra::Matrix{ComplexF64}       = zeros(ComplexF64, 0, 0)
  denoised_spectra::Matrix{Float64} = zeros(Float64, 0, 0)
  M0::Vector{Float64}               = Float64[]
  M1::Vector{Float64}               = Float64[]
  M2::Vector{Float64}               = Float64[]
  SNR::Vector{Float64}              = Float64[]
  noise_level::Vector{Float64}      = Float64[]
  a_params::Dict{Symbol, Any}       = Dict{Symbol, Any}()

  # Reference products
  ref_height::Vector{Float64}      = Float64[]
  ref_M0::Vector{Float64}          = Float64[]
  ref_M1::Vector{Float64}          = Float64[]
  ref_M2::Vector{Float64}          = Float64[]
  ref_SNR::Vector{Float64}         = Float64[]
  ref_noise_level::Vector{Float64} = Float64[]
end

Base.@kwdef mutable struct BeamSet
  vertical::RadarData = RadarData()
  north::RadarData    = RadarData()
  west::RadarData     = RadarData()
  south::RadarData    = RadarData()
  east::RadarData     = RadarData()
end

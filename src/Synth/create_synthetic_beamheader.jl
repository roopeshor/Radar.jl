function _header_overrides!(h::RadarData, overrides)
  pairs_iter = overrides isa NamedTuple ? pairs(overrides) : pairs(overrides)
  valid = Set(fieldnames(RadarData))
  for (key, value) in pairs_iter
    name = key isa Symbol ? key : Symbol(key)
    name in valid || continue
    current = getfield(h, name)
    setproperty!(h, name, convert(typeof(current), value))
  end
  return h
end

function create_synthetic_beamheader(;
  headerFields = NamedTuple(),
  altitudeRange::AbstractString = "lower",
)
  h = RadarData()
  defaults = (
    m_sMagicNumber = 369.0,
    m_sNFFT = 1024.0,
    m_sNumOfInCohIntegrations = 1.0,
    m_sCodeFlag = 0.0,
    m_sNumOfObservWindows = 1.0,
    m_sStationID = 101.0,
    m_sReserved1 = 0.0,
    m_sReserved2 = 0.0,
    m_sReserved3 = 0.0,
    m_sReserved4 = 0.0,
    m_fLatitude = 10.039999961853,
    m_fLongitude = 76.3300018310547,
    m_fAltitude = 34.0,
    m_fTotalPowerRadiated = 0.0,
    m_fRadiatedPower_ClusterWise = zeros(13),
    m_sOperationMode = 0.0,
    m_sTotalNumberofBeams = 5.0,
    m_cComment = fill(UInt8(' '), 16),
    m_fTRP_RF_Delay_us = 2.09999990463257,
    m_fMGC_Gain_dB = 30.0,
    m_usWindowType = 3.0,
    m_usReserved = 0.0,
    m_ucExperimentDataEnable = 1.0,
    m_cReserved = zeros(Int8, 3),
    m_fWindow2StartHeight = 0.0,
    m_fWindow2EndHeight = 0.0,
    m_fWindow3StartHeight = 0.0,
    m_fWindow3EndHeight = 0.0,
    m_fWindow4StartHeight = 0.0,
    m_fWindow4EndHeight = 0.0,
    m_fWindow5StartHeight = 0.0,
    m_fWindow5EndHeight = 0.0,
    m_sCurrentBeamCnt = 0.0,
  )
  _header_overrides!(h, defaults)

  if altitudeRange == "lower"
    mode = (
      173,
      0.300000011920929,
      256,
      62.5000038146973,
      0.300000011920929,
      1,
      315.0,
      8100.0,
      "0p3_30Hz",
    )
    h.m_ulCodeA = UInt32[1, 0xffffffff]
    h.m_ulCodeB = UInt32[1, 0xffffffff]
  elseif altitudeRange == "middle"
    mode = (
      94,
      1.20000004768372,
      100,
      161.290328979492,
      19.2000007629395,
      16,
      3150.0,
      20000.0,
      "1p2_30Hz",
    )
    h.m_ulCodeA = UInt32[60898, 0xffffffff]
    h.m_ulCodeB = UInt32[60701, 0xffffffff]
  else
    mode = (
      73,
      2.40000009536743,
      64,
      263.157897949219,
      38.4000015258789,
      16,
      6030.0,
      32000.0,
      "2p4_30Hz",
    )
    h.m_ulCodeA = UInt32[60898, 0xffffffff]
    h.m_ulCodeB = UInt32[60701, 0xffffffff]
  end
  h.m_sNumOfRangeBins,
  h.m_fBaudLength_us,
  h.m_sNumOfCohIntegrations,
  h.m_fIntrPulsePeriod_us,
  h.m_fPulseWidth_us,
  h.m_usCodeLength,
  h.m_fWindow1StartHeight,
  h.m_fWindow1EndHeight,
  remark = mode
  h.m_cArrRemarks =
    vcat(collect(codeunits(remark)), zeros(UInt8, 512-length(remark)))

  _header_overrides!(h, headerFields)
  h.m_fOffZenith = iszero(h.m_sCurrentBeamCnt) ? 0.0 : 10.0
  h.m_fAzimuth = get(Data.cnt2az, round(Int, h.m_sCurrentBeamCnt), 0.0)
  add_extra_RD_params!(h)
  return h
end

function get_beam_header_fields(
  b::RadarData;
  add_timestamp::Bool = false,
  remove_ctrl_char::Bool = true,
)
  kv1 = Symbol.(
    split(
      "m_sMagicNumber m_sNumOfRangeBins m_fBaudLength_us m_sNFFT m_sNumOfCohIntegrations m_sNumOfInCohIntegrations m_sCodeFlag m_fIntrPulsePeriod_us m_fPulseWidth_us m_sNumOfObservWindows m_sStationID m_sReserved1 m_sReserved2 m_sReserved3 m_sReserved4",
    ),
  )
  kts = Symbol.(split("m_sYear m_sMonth m_sDay m_sHour m_sMin m_sSec"))
  kv2 = Symbol.(
    split(
      "m_fLatitude m_fLongitude m_fAltitude m_fTotalPowerRadiated m_fRadiatedPower_ClusterWise m_sOperationMode m_sCurrentBeamCnt m_sTotalNumberofBeams m_usCodeLength m_cComment m_fTRP_RF_Delay_us m_fMGC_Gain_dB m_usWindowType m_usReserved m_cArrRemarks m_ulCodeA m_ulCodeB m_ucExperimentDataEnable m_cReserved m_fAzimuth m_fOffZenith m_fWindow1StartHeight m_fWindow1EndHeight m_fWindow2StartHeight m_fWindow2EndHeight m_fWindow3StartHeight m_fWindow3EndHeight m_fWindow4StartHeight m_fWindow4EndHeight m_fWindow5StartHeight m_fWindow5EndHeight",
    ),
  )
  keys = add_timestamp ? vcat(kv1, kts, kv2) : vcat(kv1, kv2)
  vals = Any[getfield(b, k) for k in kv1]
  if add_timestamp
    append!(vals, Any[getfield(b, k) for k in kts])
  end
  for k in kv2
    v = getfield(b, k)
    if k in (:m_cComment, :m_cArrRemarks)
      s = String(filter(!=(0x00), UInt8.(v)))
      remove_ctrl_char && (s = replace(s, r"[^\x20-\x7E]" => ""))
      push!(vals, s)
    elseif k in
           (:m_fRadiatedPower_ClusterWise, :m_ulCodeA, :m_ulCodeB, :m_cReserved)
      push!(vals, join(v, " "))
    else
      push!(vals, v)
    end
  end
  return keys, vals
end

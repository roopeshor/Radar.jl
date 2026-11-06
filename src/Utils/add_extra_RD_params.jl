function add_extra_RD_params!(inp::RadarData)
  inp.ipp_us       = inp.m_fIntrPulsePeriod_us
  inp.n_coh        = inp.m_sNumOfCohIntegrations
  inp.in_coh       = inp.m_sNumOfInCohIntegrations
  inp.start_height = inp.m_fWindow1StartHeight
  inp.end_height   = inp.m_fWindow1EndHeight
  inp.nRangeBins   = inp.m_sNumOfRangeBins
  inp.nFFT         = inp.m_sNFFT
  inp.direction    = Data.aoz2dir(az = inp.m_fAzimuth, oz = inp.m_fOffZenith)
  empty!(inp.a_params)
  return inp
end
add_extra_RD_params(inp::RadarData) = add_extra_RD_params!(deepcopy(inp))

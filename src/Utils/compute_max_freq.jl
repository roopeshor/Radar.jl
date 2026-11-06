function compute_max_freq(ipp::Real, coh::Real = 1)
  effective_sampling_time = ipp * 1e-6 * coh
  effective_sampling_time > 0 ||
    throw(ArgumentError("IPP × coherent integrations must be positive"))
  return 1 / (2.0 * effective_sampling_time)
end

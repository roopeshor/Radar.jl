"""
    compute_max_freq(ipp, coh)

Computes unambiguous Doppler frequency for a radar beam.

# Arguments
- `ipp::Real`: Inter-Pulse Period in microseconds (from Header.m_fIntrPulsePeriod_us).
- `coh::Real=1`: Number of coherent integrations (from Header.m_sNumOfCohIntegrations).

# Returns
- Maximum unambiguous Doppler frequency (Hz).
"""
function compute_max_freq(ipp::Real, coh::Real = 1)
    effective_sampling_time = ipp * 1e-6 * coh
    return 1 / (2.0 * effective_sampling_time)
end

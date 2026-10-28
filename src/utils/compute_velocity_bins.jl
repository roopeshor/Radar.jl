"""
    compute_velocity_bins(nfft, ipp, n_coh)

UTILS.compute_velocity_bins Calculates the 1D Doppler velocity axis (m/s) for a radar beam.

Evaluates Doppler frequencies fd = (-NFFT/2 : NFFT/2-1) / (NFFT * effective_PRI)
and converts to radial velocity v = - lambda * fd / 2.

# Arguments
- `nfft::Real`: BeamHeader or RadarData beam object. number of FFT points
- `ipp::Real`: m_fIntrPulsePeriod_us
- `n_coh::Real`: number of incoherent integration

# Returns
- 1 x NFFT vector of Doppler velocities (m/s).
"""
function compute_velocity_bins(nfft::Real, ipp::Real, n_coh::Real)
    c = 299792458
    fc = 205e6
    lambda = c / fc

    effective_PRI = ipp * 1e-6 * n_coh

    fd = ((-nfft/2):(nfft/2 - 1)) ./ (nfft * effective_PRI)
    vel_axis = lambda .* fd ./ 2

    return vel_axis
end

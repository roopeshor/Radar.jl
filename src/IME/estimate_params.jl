function estimate_params(beam::RadarData; vel_axis::AbstractVector)
  if isempty(beam.noise_level)
    beam.noise_level = HS_noise_estimate(beam.denoised_spectra, beam.in_coh)
  end
  candidates =
    get_candidate_peaks(beam.denoised_spectra, beam.noise_level, vel_axis)
  delta_R = 3.0e8 * beam.m_fBaudLength_us * 1e-6 / 2
  profile = trace_profile(candidates, beam.m_fOffZenith, delta_R)
  beam.M1 = profile.velocity ./ (3e8 / 205e6 / 2)
  beam.a_params[:ime_candidates] = candidates
  beam.a_params[:ime_radial_vel] = profile.velocity
  beam.a_params[:ime_power] = profile.power
  beam.a_params[:ime_snr_db] = profile.snr_db
  return beam
end

function denoise_beam(beam::RadarData)
  praw =
    isempty(beam.denoised_spectra) ? Float64.(real.(beam.spectra)) :
    beam.denoised_spectra
  nr, nfft = size(praw)
  dc = fld(nfft, 2) + 1
  if dc >= 3 && dc + 2 <= nfft
    avg = (praw[:, dc - 2] .+ praw[:, dc + 2]) ./ 2
    for j = (dc - 1):(dc + 1)
      praw[:, j] .= avg
    end
  end
  praw = moving_mean(praw, 5; dims = 2)
  noise = HS_noise_estimate(praw, max(beam.m_sNumOfInCohIntegrations, 1))
  clean = max.(praw .- reshape(noise, :, 1), 0.0)
  beam.denoised_spectra = clean
  beam.noise_level = noise
  return beam
end

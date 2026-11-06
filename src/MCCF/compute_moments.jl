function compute_moments(beam::RadarData)
  p =
    isempty(beam.denoised_spectra) ? abs2.(beam.spectra) : beam.denoised_spectra
  nr, nfft = size(p)
  m0 = zeros(Float64, nr)
  m1 = fill(NaN, nr)
  m2 = fill(NaN, nr)
  prev_idx = round(Int, nfft/2)
  costs = zeros(Float64, nr, nfft)
  for r = 1:nr
    row = @view p[r, :]
    locs = local_peaks(row)
    isempty(locs) && (locs = [argmax(row)])
    best_penalty = Inf
    best_loc = first(locs)
    for loc in locs
      peak_power = row[loc]
      ppen = 1 / (peak_power + 1e-6)
      jpen = r == 1 ? 0.0 : abs(loc-prev_idx)/nfft
      penalty = 0.8*ppen + 0.2*jpen
      costs[r, loc] = penalty
      if penalty < best_penalty
        best_penalty = penalty
        best_loc = loc
      end
    end
    prev_idx = best_loc
    left = best_loc
    right = best_loc
    while left > 1 && row[left - 1] > 0
      left -= 1
    end
    while right < nfft && row[right + 1] > 0
      right += 1
    end
    block = row[left:right]
    m0[r] = sum(block)
    m0[r] <= 0 && continue
    indices = (left:right) .- 1
    freq = (indices .- nfft/2) ./ (beam.ipp_us * beam.n_coh * nfft * 1e-6)
    m1[r] = sum(block .* freq) / m0[r]
    m2[r] = sum((freq .- m1[r]) .^ 2 .* block) / m0[r]
  end
  beam.M0 = m0
  beam.M1 = m1
  beam.M2 = m2
  beam.a_params[:mccf_cost_scores] = costs
  return beam
end

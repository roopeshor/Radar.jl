function compute_all_params(
  beam::RadarData;
  prev_frame_peaks::AbstractVector = Float64[],
  apply_mclms::Bool = true,
  win_coeff::Real = 2.0,
  wind_shear::Real = 1.0,
  velocity_res::Real = 0.06,
  cost_weights::Tuple{<:Real, <:Real} = (0.5, 0.5),
)
  data = beam.BeamData
  nr, nfft, _ = size(data)
  alpha = get(beam.a_params, :alpha, Float64(win_coeff))
  max_shear = get(beam.a_params, :max_shear, Float64(wind_shear))
  dv = get(beam.a_params, :velocity_res, Float64(velocity_res))
  df = round(Int, max_shear / dv)
  w1, w2 = cost_weights
  mclms_len = 10

  transformed = mean(fft(data, 2); dims = 3)
  spectra_complex = dropdims(transformed; dims = 3)
  noise_val = [median(abs.(@view spectra_complex[r, :])) for r = 1:nr]
  f = copy(spectra_complex)
  mask = [1, 2, 3, nfft, nfft-1]
  for r = 1:nr
    f[r, mask] .= noise_val[r] .* cis.(2π .* rand(length(mask)))
  end
  scale = max(maximum(abs, f), eps(Float64))
  cleaned = ifft(f ./ scale, 2)

  cost_scores = zeros(Float64, nr, nfft)
  denoised = zeros(Float64, nr, nfft)
  peaks_idx = zeros(Int, nr)
  st_windows = zeros(Float64, nr, 2)
  previous_peak = NaN

  for r = 1:nr
    signal = vec(cleaned[r, :])
    if r == 1
      stw = collect(1:nfft)
    else
      center = isfinite(previous_peak) ? previous_peak : nfft/2
      lo = clamp(round(Int, center-alpha*df), 1, nfft)
      hi = clamp(round(Int, center+alpha*df), 1, nfft)
      if lo > hi
        lo, hi = hi, lo
      end
      w_space = collect(lo:hi)
      if length(prev_frame_peaks) >= r && isfinite(prev_frame_peaks[r])
        t0 = prev_frame_peaks[r]
        lo2 = clamp(round(Int, t0-alpha*df), 1, nfft)
        hi2 = clamp(round(Int, t0+alpha*df), 1, nfft)
        if lo2 > hi2
          lo2, hi2 = hi2, lo2
        end
        inter = intersect(w_space, collect(lo2:hi2))
        stw = length(inter) >= alpha*df ? inter : w_space
      else
        stw = w_space
      end
    end
    st_windows[r, :] .= (first(stw), last(stw))

    seg_len = fld(nfft, 2)
    num_delays = 8
    step = fld(nfft-seg_len, num_delays-1)
    v_r1 = signal[1:seg_len]
    v_conv = Matrix{ComplexF64}(undef, num_delays-1, seg_len)
    max_powers = zeros(Float64, num_delays-1)
    for d = 2:num_delays
      start_idx = 1 + (d-1)*step
      v_ri = signal[start_idx:(start_idx + seg_len - 1)]
      v_ci = _conv1d(v_r1, conj.(reverse(v_ri)))
      lo = fld(seg_len, 2)
      hi = fld(3seg_len, 2)-1
      center = v_ci[lo:hi]
      v_conv[d - 1, :] .= center
      spec_d = abs2.(fftshift(fft_n(center, nfft)))
      max_powers[d - 1] = maximum(spec_d[stw])
    end
    best_d = argmax(max_powers)
    x_opt = vec(v_conv[best_d, :])
    adjacent = best_d > 1 ? best_d-1 : min(2, num_delays-1)
    x_sub = vec(v_conv[adjacent, :])
    y = apply_mclms ? mclms_filter(x_opt, x_sub, mclms_len) : x_opt
    spec_mclms = abs2.(fftshift(fft_n(y, nfft)))
    sm = maximum(spec_mclms)
    sm > 0 && (spec_mclms ./= sm)
    denoised[r, :] .= spec_mclms

    if r == 1
      avg_conv = vec(mean(v_conv; dims = 1))
      final_spec = abs2.(fftshift(fft_n(avg_conv, nfft)))
      selected = argmax(final_spec)
    else
      window_spec = spec_mclms[stw]
      max_power = maximum(window_spec)
      peak_locs = local_peaks(window_spec)
      valid_locs = filter(i -> window_spec[i] >= 0.4max_power, peak_locs)
      if isempty(valid_locs)
        valid_locs = [argmax(window_spec)]
      end
      order = sortperm(window_spec[valid_locs]; rev = true)
      valid_locs = valid_locs[order[1:min(3, length(order))]]
      candidate_indices = stw[valid_locs]
      candidate_powers = window_spec[valid_locs]
      best_cost = -Inf
      selected = first(candidate_indices)
      for (idx, pwr) in zip(candidate_indices, candidate_powers)
        vcost = 1 - abs((idx-previous_peak)*dv) / (alpha*max_shear)
        cost = w1*vcost + w2*pwr
        cost_scores[r, idx] = cost
        if cost > best_cost
          best_cost = cost
          selected = idx
        end
      end
    end
    previous_peak = selected
    peaks_idx[r] = selected
    left = selected
    right = selected
    while left > 1 && spec_mclms[left - 1] > 0
      left -= 1
    end
    while right < nfft && spec_mclms[right + 1] > 0
      right += 1
    end
  end

  fmax = (nfft/2) / (nfft * beam.ipp_us * 1e-6 * beam.n_coh)
  beam.M1 = vec(rescale_linear(peaks_idx, -fmax, fmax, 1, nfft))
  beam.denoised_spectra = denoised
  beam.a_params[:mccf_cost_scores] = cost_scores
  beam.a_params[:mccf_ST_W] = st_windows
  beam.a_params[:peaks_idx] = peaks_idx
  return beam
end

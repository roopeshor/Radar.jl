struct PeakCandidate
  idx::Int
  velocity::Float64
  power::Float64
  snr_db::Float64
  width::Int
end

function get_candidate_peaks(
  pclean::AbstractMatrix,
  noise_level::AbstractVector,
  vel_axis::AbstractVector,
)
  nr, nfft = size(pclean)
  length(noise_level) == nr || throw(
    DimensionMismatch("noise vector length must equal number of range bins"),
  )
  length(vel_axis) == nfft ||
    throw(DimensionMismatch("velocity axis length must equal NFFT"))
  candidates = [PeakCandidate[] for _ = 1:nr]
  half_width = 0.10 * (maximum(vel_axis) - minimum(vel_axis))
  previous_velocity = NaN

  for r = 1:nr
    p = Float64.(pclean[r, :])
    sn = noise_level[r]
    center = if r == 1
      vel_axis[argmax(p)]
    else
      previous_velocity
    end
    if !isfinite(center)
      center = vel_axis[argmax(p)]
    end
    idx_window = findall(
      (vel_axis .>= center-half_width) .& (vel_axis .<= center+half_width),
    )
    isempty(idx_window) && continue
    working = copy(p[idx_window])
    selected = PeakCandidate[]
    for _ = 1:3
      isempty(working) && break
      max_amp, local_idx = findmax(working)
      max_amp <= 0 && break
      left = local_idx
      while left > 1 && working[left - 1] > 0
        left -= 1
      end
      right = local_idx
      while right < length(working) && working[right + 1] > 0
        right += 1
      end
      ppeak = working[left:right]
      vpeak = vel_axis[idx_window[left:right]]
      total_power = sum(ppeak)
      total_power <= 0 && break
      mean_velocity = sum(ppeak .* vpeak) / total_power
      snr_linear = sn > 0 ? total_power / (sn*nfft) : Inf
      snr_db = 10log10(snr_linear)
      snr_db <= 0 && break
      push!(
        selected,
        PeakCandidate(
          idx_window[local_idx],
          mean_velocity,
          total_power,
          snr_db,
          right-left+1,
        ),
      )
      working[left:right] .= 0
    end
    candidates[r] = selected
    if !isempty(selected)
      closest = argmin(abs.([c.velocity for c in selected] .- center))
      previous_velocity = selected[closest].velocity
    elseif r == 1
      previous_velocity = center
    end
  end
  return candidates
end

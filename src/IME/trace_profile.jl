function trace_profile(
  candidates::Vector{Vector{PeakCandidate}},
  off_zenith::Real,
  delta_R::Real,
)
  nr = length(candidates)
  velocity = zeros(Float64, nr)
  power = zeros(Float64, nr)
  snr_db = zeros(Float64, nr)
  w1 = 0.4
  w2 = 0.6
  group_size = 5
  delta_vm = max(0.5, 0.05 * sind(off_zenith) * delta_R)

  for start_bin = 1:group_size:nr
    end_bin = min(start_bin + group_size - 1, nr)
    bins = end_bin-start_bin+1
    peaks0 = candidates[start_bin]
    if isempty(peaks0)
      candidates[start_bin] = [PeakCandidate(0, NaN, NaN, NaN, 0)]
      peaks0 = candidates[start_bin]
    end
    costs = Vector{Vector{Float64}}(undef, bins)
    paths = Vector{Vector{Vector{Int}}}(undef, bins)
    snrmax = maximum(c.snr_db for c in peaks0)
    snrmax <= 0 && (snrmax = 1.0)
    costs[1] = [w2 * c.snr_db / snrmax for c in peaks0]
    paths[1] = [[p] for p in eachindex(peaks0)]

    for localbin = 2:bins
      r = start_bin + localbin - 1
      prev = candidates[r - 1]
      curr = candidates[r]
      isempty(curr) && (
        candidates[r] = [PeakCandidate(0, NaN, NaN, NaN, 0)];
        curr = candidates[r]
      )
      snrmax = maximum(c.snr_db for c in curr)
      snrmax <= 0 && (snrmax = 1.0)
      newcost = fill(-Inf, length(curr))
      newpaths = Vector{Vector{Int}}(undef, length(curr))
      for cidx in eachindex(curr)
        c = curr[cidx]
        best_cost = -Inf
        best_prev = 1
        pterm = c.snr_db / snrmax
        for pidx in eachindex(prev)
          p = prev[pidx]
          x =
            (!isfinite(c.velocity) || !isfinite(p.velocity)) ? Inf :
            abs(c.velocity-p.velocity)
          tx =
            x <= delta_vm ? 1.0 :
            (x <= 2delta_vm ? (2delta_vm-x)/delta_vm : 0.0)
          total = costs[localbin - 1][pidx] + w1*tx + w2*pterm
          if total > best_cost
            best_cost = total
            best_prev = pidx
          end
        end
        newcost[cidx] = best_cost
        newpaths[cidx] = [paths[localbin - 1][best_prev]; cidx]
      end
      costs[localbin] = newcost
      paths[localbin] = newpaths
    end
    best = argmax(costs[bins])
    path = paths[bins][best]
    for localbin = 1:bins
      r = start_bin + localbin - 1
      c = candidates[r][path[localbin]]
      velocity[r] = c.velocity
      power[r] = c.power
      snr_db[r] = c.snr_db
    end
  end
  return (; velocity, power, snr_db)
end

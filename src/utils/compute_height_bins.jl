"""
    compute_height_bins(start_range, end_range, num_range_bins)

UTILS.compute_height_bins Computes linearly spaced altitude bin values across observation windows.

# Arguments
- `start_range::Real`: Starting observation window height (km or m).
- `end_range::Real`: Ending observation window height (km or m).
- `num_range_bins::Real`: Total number of range bins.

# Returns
- Vector of altitude range bin heights.
"""
function compute_height_bins(start_range::Real, end_range::Real, num_range_bins::Real)
    n = Int(num_range_bins)
    return range(start_range, end_range, length=n) |> collect
end

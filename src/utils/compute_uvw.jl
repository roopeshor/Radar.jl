"""
    compute_uvw(obs)

Computes uvw from moments.

# Arguments
- `obs::Observation`

# Returns
- `obs::Observation`
"""
function compute_uvw(obs::Observation)
    sum_M1 = obs.east.M1 .+ obs.west.M1 .+ obs.north.M1 .+ obs.south.M1
    cth = cosd(obs.north.m_fOffZenith)

    obs.U = .-(obs.east.M1 .- obs.west.M1) .* obs.DBS_Factor_H
    obs.V = .-(obs.north.M1 .- obs.south.M1) .* obs.DBS_Factor_H
    obs.W = .-obs.DBS_Factor .* (cth .* sum_M1 .+ obs.vertical.M1) ./ (4 * cth^2 + 1)

    return obs
end

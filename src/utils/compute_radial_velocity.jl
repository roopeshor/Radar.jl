"""
    compute_radial_velocity(u, v, w, az, oz)

Projects 3D wind velocity vector (u, v, w) along the radar beam direction.

Calculates line-of-sight radial velocity vr = (u*sin(az) + v*cos(az))*sin(oz) + w*cos(oz)
given azimuth angle az (degrees) and off-zenith angle oz (degrees).

# Arguments
- `u::AbstractVector{<:Real}`: Zonal wind component profile (m/s).
- `v::AbstractVector{<:Real}`: Meridional wind component profile (m/s).
- `w::AbstractVector{<:Real}`: Vertical wind component profile (m/s).
- `az::Real`: Beam azimuth angle in degrees (0=North, 90=East, 180=South, 270=West).
- `oz::Real`: Beam off-zenith angle in degrees.

# Returns
- Projected line-of-sight radial velocity profile (m/s).
"""
function compute_radial_velocity(u::AbstractVector{<:Real}, v::AbstractVector{<:Real}, w::AbstractVector{<:Real}, az::Real, oz::Real)
    if oz == 0
        return w
    else
        return (u .* sind(az) .+ v .* cosd(az)) .* sind(oz) .+ w .* cosd(oz)
    end
end

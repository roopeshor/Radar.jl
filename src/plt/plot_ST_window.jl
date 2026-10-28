"""
    plot_ST_window(stw, v_max, heights, box_color="#f0f")

# Arguments
- `stw::AbstractMatrix{<:Real}`: ST window, array of tuple of (minIDX, maxIDX)
- `v_max::Real`: max velocity
- `heights::AbstractVector{<:Real}`: height bins
- `box_color::Any="#f0f"`
"""
function plot_ST_window(
    stw::AbstractMatrix{<:Real},
	v_max::Real,
	heights::AbstractVector{<:Real},
	box_color::Any="#f0f"
)
    stw_rescaled = -v_max .+ (stw .- 1.0) .* (2.0 * v_max) ./ (1024.0 - 1.0)
    h = heights[2] - heights[1]

    for i in 1:length(heights)
        x = stw_rescaled[i, 1]
        y = heights[i]
        w = stw_rescaled[i, 2] - x

        plot!([-x - w, -x], [y, y], linewidth=1, color=box_color, label=false)
        plot!([-x - w, -x - w], [y - h / 3, y + h / 3], linewidth=1, color=box_color, label=false)
        plot!([-x, -x], [y - h / 3, y + h / 3], linewidth=1, color=box_color, label=false)
    end
end

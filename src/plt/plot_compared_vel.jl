"""
    plot_compared_vel(; calc, ref, plotidx, heights, title="", subplots=3)

plt.plot_compared_vel Plots calculated vs reference wind velocity profiles against altitude.

# Arguments
- `calc::AbstractVector{<:Real}`: Computed velocity profile vector.
- `ref::AbstractVector{<:Real}`: Reference velocity profile vector.
- `plotidx::Integer`: Subplot index.
- `heights::AbstractVector{<:Real}`: Altitude heights array.
- `title::AbstractString=""`: Subplot title string.
- `subplots::Integer=3`: Total subplots count.
"""
function plot_compared_vel(; calc::AbstractVector{<:Real},
    ref::AbstractVector{<:Real},
    plotidx::Integer,
    heights::AbstractVector{<:Real},
    title::AbstractString="",
    subplots::Integer=3
)
    plot(ref, heights, label="ref", layout=(1, subplots), subplot=plotidx)
    plot!(calc, heights, label="calc", subplot=plotidx)

    xlims!(minimum(ref) * 2, maximum(ref) * 2)
    ylims!(minimum(heights), maximum(heights))
    xlabel!("wind velocity (m/s)")
    ylabel!("height (km)")
    title!(title)
end

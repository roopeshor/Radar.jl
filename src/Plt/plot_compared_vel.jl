function plot_compared_vel(;
  calc::AbstractVector,
  ref::AbstractVector,
  plotidx::Integer = 1,
  heights::AbstractVector,
  title::AbstractString = "",
  subplots::Integer = 3,
  plot_handle = Plots.current(),
)
  p = plot!(plot_handle, ref, heights; label = "reference", subplot = plotidx)
  plot!(p, calc, heights; label = "computed", subplot = plotidx)
  xlabel!(p, "wind velocity (m/s)")
  ylabel!(p, "height (km)")
  title!(p, title)
  if !isempty(ref)
    xlims!(p, (minimum(ref)*2, maximum(ref)*2))
  end
  !isempty(heights) && ylims!(p, (minimum(heights), maximum(heights)))
  return p
end

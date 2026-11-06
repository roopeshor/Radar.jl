function plot_ST_window(
  stw::AbstractMatrix,
  v_max::Real,
  heights::AbstractVector,
  box_color = "magenta";
  plot_handle = Plots.current(),
)
  size(stw, 2) == 2 ||
    throw(DimensionMismatch("ST window must have two columns"))
  length(heights) == size(stw, 1) ||
    throw(DimensionMismatch("height count must match ST windows"))
  mapped = (stw .- 1) .* (2v_max/1023) .- v_max
  dh = length(heights) > 1 ? heights[2]-heights[1] : 1.0
  for i in eachindex(heights)
    x1 = -mapped[i, 2]
    x2 = -mapped[i, 1]
    y = heights[i]
    plot!(
      plot_handle,
      [x1, x2],
      [y, y];
      color = box_color,
      linewidth = 1,
      label = false,
    )
    plot!(
      plot_handle,
      [x1, x1],
      [y-dh/3, y+dh/3];
      color = box_color,
      linewidth = 1,
      label = false,
    )
    plot!(
      plot_handle,
      [x2, x2],
      [y-dh/3, y+dh/3];
      color = box_color,
      linewidth = 1,
      label = false,
    )
  end
  return plot_handle
end

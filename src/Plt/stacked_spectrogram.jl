function stacked_spectrogram(
  spectrum::AbstractMatrix,
  y_ticks::AbstractVector = collect(1:size(spectrum, 1)),
  x_ticks::AbstractVector = collect(1:size(spectrum, 2));
  amplitude_scale::Real = 1.2,
  plot_color = "blue",
  plot_gradiated::Bool = false,
  plot_add_colorbar::Bool = false,
  plot_thick::Real = 1.0,
  baseline_color = "gray",
  baseline_thick::Real = 0.0,
  peak_marker = "blue*",
  marker_size::Real = 4,
  ref_line_pos::Real = 0,
  ref_line_thick::Real = 1,
  ref_line_color = "red",
  decimate_factor::Integer = 1,
  x_axis_type::AbstractString = "velocity",
  plot_handle = nothing,
)
  decimate_factor > 1 && (
    x_ticks = x_ticks[1:decimate_factor:end];
    spectrum = spectrum[:, 1:decimate_factor:end]
  )
  mat =
    x_axis_type == "velocity" ? reverse(Float64.(spectrum); dims = 2) :
    Float64.(spectrum)
  dy = length(y_ticks) > 1 ? y_ticks[2]-y_ticks[1] : 1.0
  finitevals = filter(isfinite, vec(mat))
  ydiff = isempty(finitevals) ? 1.0 : maximum(finitevals)-minimum(finitevals)
  scale = ydiff == 0 ? amplitude_scale*dy : amplitude_scale*dy/ydiff
  p = plot_handle === nothing ? plot() : plot_handle
  if ref_line_thick > 0
    vline!(
      p,
      [ref_line_pos];
      color = ref_line_color,
      linewidth = ref_line_thick,
      label = false,
    )
  end
  for r in axes(mat, 1)
    yy = y_ticks[r] .+ mat[r, :] .* scale
    col = plot_gradiated ? cgrad(:viridis)[r / size(mat, 1)] : plot_color
    plot!(p, x_ticks, yy; color = col, linewidth = plot_thick, label = false)
    baseline_thick > 0 && hline!(
      p,
      [y_ticks[r]];
      color = baseline_color,
      linewidth = baseline_thick,
      label = false,
    )
    if !isempty(peak_marker) && marker_size > 0
      j = argmax(mat[r, :])
      scatter!(
        p,
        [x_ticks[j]],
        [y_ticks[r]+mat[r, j]*scale];
        marker = :circle,
        markersize = marker_size,
        color = col,
        label = false,
      )
    end
  end
  xlabel!(p, "Doppler $x_axis_type (m/s)")
  ylabel!(p, "Altitude (km)")
  xlims!(p, (minimum(x_ticks), maximum(x_ticks)))
  ylims!(p, (minimum(y_ticks)-scale, maximum(y_ticks)+scale))
  return p
end

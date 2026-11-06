function plot_doppler_spectra(;
  spectra::AbstractMatrix,
  comp_m::AbstractVector = Float64[],
  ref_m::AbstractVector = Float64[],
  direction::AbstractString = "",
  heights::AbstractVector = collect(1:size(spectra, 1)),
  x_max::Real = 30.0,
  add_legend::Bool = true,
  comp_m_color = "black",
  ref_m_color = "red",
  plot_thickness::Real = 1.0,
  x_axis_type::AbstractString = "velocity",
  plot_handle = nothing,
)
  mat = x_axis_type == "velocity" ? reverse(spectra; dims = 2) : spectra
  x = collect(range(-x_max, x_max, length = size(mat, 2)))
  params = Dict(
    :yflip    => false,
    :xlabel   => "Doppler $x_axis_type (m/s)",
    :ylabel   => "Altitude (km)",
    :title    => direction,
    :colorbar => true,
    :legend   => add_legend,
    :size     => (450, 600),
    :dpi      => 200,
    :c        => PARULA,
  )
  p = if plot_handle === nothing
    heatmap(x, heights, mat; params...)
  else
    heatmap!(plot_handle, x, heights, mat; params...)
  end
  if !isempty(ref_m)
    plot!(
      p,
      ref_m,
      heights;
      color = ref_m_color,
      linewidth = plot_thickness,
      label = "reference",
    )
  end
  if !isempty(comp_m)
    plot!(
      p,
      comp_m,
      heights;
      color = comp_m_color,
      linewidth = plot_thickness,
      label = "computed",
    )
  end
  xlims!(p, (-x_max, x_max))
  return p
end

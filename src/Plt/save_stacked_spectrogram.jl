"Render and save the stacked spectrogram using the Julia plotting backend."
function save_stacked_spectrogram(
  spectrum::AbstractMatrix,
  y_ticks::AbstractVector = collect(1:size(spectrum, 1)),
  x_ticks::AbstractVector = collect(1:size(spectrum, 2));
  export_path::AbstractString = "plots/stacked_spectrogram.png",
  export_dpi::Real = 300,
  kwargs...,
)
  mkpath(dirname(export_path))
  p = stacked_spectrogram(spectrum, y_ticks, x_ticks; kwargs...)
  plot!(p; dpi = export_dpi)
  savefig(p, export_path)
  return export_path
end

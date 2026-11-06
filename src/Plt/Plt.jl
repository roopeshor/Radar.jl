module Plt
using Colors

# if doing benchmarking, dont show plot with:
# ENV["GKSwstype"] = "100"
using Plots

export plot_ST_window,
  plot_compared_vel,
  plot_doppler_spectra,
  save_stacked_spectrogram,
  stacked_spectrogram

include("colormaps.jl")
include("plot_ST_window.jl")
include("plot_compared_vel.jl")
include("plot_doppler_spectra.jl")
include("stacked_spectrogram.jl")
include("save_stacked_spectrogram.jl")

end # module Plt

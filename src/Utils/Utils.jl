module Utils

using Statistics
using Random
using FFTW
using ..Types: RadarData, Observation, Data

export add_extra_RD_params!,
  add_extra_RD_params,
  compute_height_bins,
  compute_max_freq,
  compute_radial_velocity,
  compute_uvw,
  compute_velocity_bins,
  get_beam_header_fields,
  get_original_spectrum,
  fftshift,
  ifftshift,
  hanning_window,
  moving_mean,
  gaussian_smooth,
  local_peaks,
  rescale_linear,
  safe_corrcoef,
  fft_n,
	gradient

include("helpers.jl")
include("add_extra_RD_params.jl")
include("compute_height_bins.jl")
include("compute_max_freq.jl")
include("compute_radial_velocity.jl")
include("compute_uvw.jl")
include("compute_velocity_bins.jl")
include("get_beam_header_fields.jl")
include("get_original_spectrum.jl")

end # module Utils

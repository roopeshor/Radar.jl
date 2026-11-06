module MCCF

using Statistics
using Random
using FFTW
using ..Types: RadarData
using ..Utils: fftshift, gaussian_smooth, local_peaks, rescale_linear, fft_n

export compute_all_params,
  compute_beam_spectra,
  compute_beam_spectra_vec,
  compute_mc_spectra,
  compute_moments,
  mclms_filter

include("helpers.jl")
include("mclms_filter.jl")
include("compute_beam_spectra.jl")
include("compute_beam_spectra_vec.jl")
include("compute_mc_spectra.jl")
include("compute_moments.jl")
include("compute_all_params.jl")

end # module MCCF

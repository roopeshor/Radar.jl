module IME

using Statistics
using FFTW
using ..Types: RadarData
using ..Utils: hanning_window, fftshift, moving_mean

export HS_noise_estimate,
  calculate_moments,
  compute_beam_spectra,
  denoise_beam,
  estimate_params,
  get_candidate_peaks,
  trace_profile

include("HS_noise_estimate.jl")
include("calculate_moments.jl")
include("compute_beam_spectra.jl")
include("denoise_beam.jl")
include("get_candidate_peaks.jl")
include("trace_profile.jl")
include("estimate_params.jl")

end # module IME

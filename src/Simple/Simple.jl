module Simple

using Statistics
using FFTW
using ..Types: RadarData
using ..Utils: fftshift, moving_mean

export HS_noise_estimate, compute_beam_spectra, compute_moments, denoise_beam

include("HS_noise_estimate.jl")
include("compute_beam_spectra.jl")
include("compute_moments.jl")
include("denoise_beam.jl")

end # module Simple

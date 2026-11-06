module Synth

using Statistics
using Random
using FFTW
using ..Types: RadarData, Observation, Data
import ..Utils
using ..Utils:
  gaussian_smooth,
  add_extra_RD_params!,
  compute_radial_velocity,
  compute_height_bins,
  compute_velocity_bins,
  compute_max_freq,
  get_original_spectrum
using ..RadarIO: write_raw_file

@inline _round_away(x::Real) = copysign(floor(abs(x) + 0.5), x)

export compute_spectrum_broadening,
  convective_plumes,
  create_synthetic_beamheader,
  create_synthetic_observation,
  generate_wind_profile,
  gravity_waves,
  kolmogorov_turbulence,
  orographic_updraft,
  somali_jet,
  stratospheric_qbo,
  synoptic_bg,
  synthesize_iq_data,
  tropical_easterly_jet

include("compute_spectrum_broadening.jl")
include("convective_plumes.jl")
include("create_synthetic_beamheader.jl")
include("gravity_waves.jl")
include("kolmogorov_turbulence.jl")
include("orographic_updraft.jl")
include("somali_jet.jl")
include("stratospheric_qbo.jl")
include("synoptic_bg.jl")
include("tropical_easterly_jet.jl")
include("generate_wind_profile.jl")
include("synthesize_iq_data.jl")
include("create_synthetic_observation.jl")

end # module Synth

# Radar.jl - Top-level module for the ST Radar signal processing toolbox.
#
# Load order:
#   1. Data, DynamicProps, BeamHeader, RadarData  (shared types)
#   2. Submodules that only need RadarData (io, utils, simple, mccf, ime, plt)
#   3. Observation.jl  (defines Observation struct + apply_beam_fx + compute_uvw)
#   4. synth  (needs Observation, so comes last)

module Radar

include("Data.jl")
include("DynamicProps.jl")
include("BeamHeader.jl")
include("RadarData.jl")

include("io.jl")
include("Observation.jl")

# ── io ─────────────────────────────────────────────────────────────────────
# module io
#     using ..Radar: RadarData, BeamHeader, DynamicProps, Data, alias_params!
#     using CSV, DataFrames
# end

# ── utils ──────────────────────────────────────────────────────────────────
# module utils
#     using ..Radar: RadarData, BeamHeader
#     using Statistics: mean
#     using FFTW: fft, fftshift
#     include("utils/compute_height_bins.jl")
#     include("utils/compute_max_freq.jl")
#     include("utils/compute_radial_velocity.jl")
#     include("utils/compute_velocity_bins.jl")
#     include("utils/get_beam_header_fields.jl")
#     include("utils/get_original_spectrum.jl")
#     # compute_uvw is defined in Observation.jl (needs Observation type)
# end

# # ── simple ─────────────────────────────────────────────────────────────────
# module simple
#     using ..Radar: RadarData, BeamHeader
#     using Statistics: mean, median
#     using FFTW: fft, fftshift
#     include("simple/HS_noise_estimate.jl")
#     include("simple/compute_beam_spectra.jl")
#     include("simple/denoise_beam.jl")
#     include("simple/compute_moments.jl")
# end

# # ── mccf ───────────────────────────────────────────────────────────────────
# module mccf
#     using ..Radar: RadarData, BeamHeader
#     using Statistics: mean, median
#     using FFTW: fft, ifft, fftshift
#     using DSP: conv
#     include("mccf/mclms_filter.jl")
#     include("mccf/compute_mc_spectra.jl")
#     include("mccf/compute_beam_spectra.jl")
#     include("mccf/compute_beam_spectra_vec.jl")
#     include("mccf/compute_moments.jl")
#     include("mccf/compute_all_params.jl")
# end

# # ── ime ────────────────────────────────────────────────────────────────────
# module ime
#     using ..Radar: RadarData, BeamHeader
#     using Statistics: mean
#     using FFTW: fft, fftshift
#     include("ime/HS_noise_estimate.jl")
#     include("ime/compute_beam_spectra.jl")
#     include("ime/denoise_beam.jl")
#     include("ime/get_candidate_peaks.jl")
#     include("ime/calculate_moments.jl")
#     include("ime/trace_profile.jl")
#     include("ime/estimate_params.jl")
# end

# # ── plt ────────────────────────────────────────────────────────────────────
# module plt
#     using Plots
#     include("plt/plot_doppler_spectra.jl")
#     include("plt/plot_compared_vel.jl")
#     include("plt/stacked_spectrogram.jl")
#     include("plt/plot_ST_window.jl")
# end

# # ── Observation (defines struct + apply_beam_fx + compute_uvw) ─────────────
# # Included at module scope so Observation is a top-level Radar type.
# include("Observation.jl")

# # ── synth ──────────────────────────────────────────────────────────────────
# # Loaded AFTER Observation.jl so the Observation type is available in scope.
# module synth
#     import ..Radar: RadarData, BeamHeader, Data, Observation, alias_params!
#     import ..Radar.io: write_raw_file
#     import ..Radar.utils: compute_radial_velocity
#     using Statistics: mean
#     using FFTW: fft, ifft, fftshift, ifftshift
#     include("synth/compute_spectrum_broadening.jl")
#     include("synth/generate_wind_profile.jl")
#     include("synth/synthesize_iq_data_vectrzd.jl")
#     include("synth/create_synthetic_beamheader.jl")
#     include("synth/create_synthetic_observation.jl")
# end

# # Re-export commonly used symbols
# export Observation, RadarData, BeamHeader, DynamicProps, Data
# export apply_beam_fx, compute_uvw

end

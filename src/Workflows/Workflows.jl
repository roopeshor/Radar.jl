module Workflows

using Plots
using ..Types: Observation, RadarData, Data
import ..Utils
import ..RadarIO
import ..Synth

export get_observation, getObservation, apply_beam_fx, get_plotted_wind_profile

include("get_observation.jl")
include("apply_beam_fx.jl")
include("get_plotted_wind_profile.jl")

end # module Workflows

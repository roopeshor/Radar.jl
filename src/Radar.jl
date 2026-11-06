module Radar

include("Types/Types.jl")
include("Utils/Utils.jl")
include("RadarIO/RadarIO.jl")
include("IME/IME.jl")
include("MCCF/MCCF.jl")
include("Simple/Simple.jl")
include("Plt/Plt.jl")
include("Synth/Synth.jl")
include("Workflows/Workflows.jl")

const Data = Types.Data
const RadarData = Types.RadarData
const Observation = Types.Observation

export Types, Data, RadarData, Observation
export Utils, RadarIO, IME, MCCF, Simple, Plt, Synth, Workflows
export get_observation, getObservation, apply_beam_fx, get_plotted_wind_profile

using .Workflows:
  get_observation, getObservation, apply_beam_fx, get_plotted_wind_profile

end

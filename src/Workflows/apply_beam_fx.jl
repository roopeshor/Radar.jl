"Apply a beam-processing function independently to all five beams in an observation."
function apply_beam_fx(obs::Observation, fx::Function; kwargs...)
  for direction in Data.directions
    beam = getproperty(obs, direction)
    result = fx(beam; kwargs...)
    result isa RadarData ||
      throw(ArgumentError("Beam function $(fx) must return RadarData"))
    setproperty!(obs, Symbol(direction), result)
  end
  return obs
end

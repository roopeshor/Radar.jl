"Mutable collection of the five radar beams and observation-level metadata."
Base.@kwdef mutable struct Observation
  north::RadarData    = RadarData()
  east::RadarData     = RadarData()
  west::RadarData     = RadarData()
  south::RadarData    = RadarData()
  vertical::RadarData = RadarData()

  raw_filepath::String = ""

  ref_height::Vector{Float64} = Float64[]
  ref_U::Vector{Float64}      = Float64[]
  ref_V::Vector{Float64}      = Float64[]
  ref_W::Vector{Float64}      = Float64[]

  U::Vector{Float64} = Float64[]
  V::Vector{Float64} = Float64[]
  W::Vector{Float64} = Float64[]

  heights::Vector{Float64}    = Float64[]
  velocities::Vector{Float64} = Float64[]

  DBS_Factor_H::Float64 = NaN
  DBS_Factor::Float64   = NaN
  f_max::Float64        = NaN
  v_max::Float64        = NaN
  start_height::Float64 = NaN
  end_height::Float64   = NaN
  nRangeBins::Float64   = NaN
  ipp_us::Float64       = NaN
  in_coh::Float64       = NaN
  n_coh::Float64        = NaN
  nFFT::Float64         = NaN
end

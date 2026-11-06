using Radar
using Statistics
using Plots

function run_simulated_data(;
  save_plots = true,
  output_dir = joinpath("outputs", "resim", "plots"),
)
  mkpath(output_dir)
  altitude = "lower"
  headers = [
    Synth.create_synthetic_beamheader(
      headerFields = (m_sCurrentBeamCnt = i-1,),
      altitudeRange = altitude,
    ) for i = 1:5
  ]
  h, u, v, w = Synth.generate_wind_profile(
    nRangeBins = Int(headers[1].nRangeBins),
    start_height = 0,
    end_height = 20,
  )
  w .= 0
  bands = [
    8 5;
    9 5;
    10 9.4;
    11 5;
    50 5;
    51 6.2;
    52 7.4;
    53 14.2;
    95 5;
    96 5;
    97 5;
    98 5;
    134 5;
    135 5;
    136 5;
    137 17.5
  ]
  obs = Synth.create_synthetic_observation(
    u,
    v,
    w;
    headerFields = headers,
    SNR          = 10 .^ (collect(range(50, -35, length = length(h))) ./ 10),
    dc_rel_power = 0.2,
    gc_rel_power = 0.003687,
    bright_bands = Dict(d=>bands for d in ("north", "south", "east", "west")),
  )
  apply_beam_fx(obs, MCCF.compute_all_params; wind_shear = 1, win_coeff = 2)
  for direction in Data.directions_str
    beam = getproperty(obs, Symbol(direction))
    if save_plots
      p = Plt.plot_doppler_spectra(
        spectra   = 20 .* log10.(abs.(beam.spectra)),
        direction = direction,
        ref_m     = (-beam.ref_M1 .* obs.DBS_Factor),
        comp_m    = (-beam.M1 .* obs.DBS_Factor),
        heights   = obs.heights ./ 1000,
        x_max     = obs.v_max,
      )
      savefig(p, joinpath(output_dir, "1000_mccf_$(direction).png"))
    end
  end
  return obs
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
  run_simulated_data()
end

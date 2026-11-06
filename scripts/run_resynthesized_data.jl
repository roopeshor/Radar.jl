using Radar
using Plots

function run_resynthesized_data(filepath = joinpath(
  "Data",
  "other",
  "EXP_DBS_CH4_29Jul2026_16_23_15.raw",
))
  obs = get_observation(filepath)
  apply_beam_fx(obs, Simple.compute_beam_spectra)
  headers = [obs.north, obs.south, obs.east, obs.west, obs.vertical]
  synthetic = Synth.create_synthetic_observation(
    obs.ref_U,
    obs.ref_V,
    obs.ref_W;
    headerFields = headers,
  )
  apply_beam_fx(synthetic, MCCF.compute_beam_spectra)
  for direction in Data.directions
    raw = getproperty(obs, direction)
    syn = getproperty(synthetic, direction)
    p1 = Plt.plot_doppler_spectra(
      spectra = log10.(abs.(raw.spectra)),
      ref_m = (-raw.ref_M1 .* obs.DBS_Factor),
      direction = raw.direction,
      heights = obs.heights,
      x_max = obs.v_max,
    )
    p2 = Plt.plot_doppler_spectra(
      spectra = log10.(abs.(syn.spectra)),
      ref_m = (-syn.ref_M1 .* obs.DBS_Factor),
      direction = syn.direction,
      heights = synthetic.heights,
      x_max = synthetic.v_max,
    )
    display(plot(p1, p2; layout = (1, 2), title = direction))
  end
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
  run_resynthesized_data()
end

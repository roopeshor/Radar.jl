using Radar

function run_plot_spectrograph(filepath = joinpath(
  "Data",
  "rain",
  "EXP_DBS_CH4_28Jun2026_14_28_16",
))
  obs = get_observation(filepath; add_mmts = false, add_uvw = false)
  apply_beam_fx(obs, MCCF.compute_beam_spectra)
  apply_beam_fx(obs, Simple.denoise_beam)
  data = obs.north.denoised_spectra
  Plt.plot_doppler_spectra(
    spectra = data,
    heights = obs.heights,
    x_axis_type = "doppler",
    x_max = maximum(abs, obs.velocities),
  ) |> display
  Plt.save_stacked_spectrogram(
    data,
    obs.heights,
    obs.velocities;
    plot_gradiated = true,
    plot_add_colorbar = true,
    plot_color = :viridis,
    plot_thick = 0.5,
    peak_marker = "white.",
    marker_size = 0,
    amplitude_scale = 0.9,
    ref_line_thick = 0,
    ref_line_color = "red",
    export_path = joinpath("plots", "28Jun2026_14_28_16-north-raw.png"),
    export_dpi = 700,
  )
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
  run_plot_spectrograph()
end

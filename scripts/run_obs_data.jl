using Radar
using Plots

function main(filepath = joinpath(
  "Data",
  "Data_SOE",
  "Radar_data",
  "20_Jan_2026",
  "EXP_DBS_CH4_20Jan2026_14_22_11.raw",
))
  println("Processing file: ", filepath)
  obs = get_observation(filepath)
  apply_beam_fx(obs, MCCF.compute_all_params)
  Utils.compute_uvw(obs)
  d = obs.south
  p = Plt.plot_doppler_spectra(
    spectra = 10 .* log10.(abs.(d.spectra) .^ 2),
    direction = d.direction,
    ref_m = (-d.ref_M1 .* obs.DBS_Factor),
    comp_m = (-d.M1 .* obs.DBS_Factor),
    heights = obs.heights ./ 1000,
    x_max = obs.v_max,
  )
  display(p)
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
  main()

  println("\nPress Enter to exit...")
  readline()
end

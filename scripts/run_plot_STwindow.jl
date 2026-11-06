using Radar
using Plots

function run_plot_STwindow(filepath = joinpath(
  "Data",
  "other",
  "EXP_DBS_CH4_29Jul2026_16_16_17",
))
  obs = get_observation(filepath)
  apply_beam_fx(obs, MCCF.compute_all_params; win_coeff = 1.5)
  beam = obs.north
  vmax = obs.v_max
  vl = [-vmax, vmax] ./ 3
  vidx = round.(Int, Utils.rescale_linear(vl, 1, obs.nFFT, -vmax, vmax))
  vidx = clamp.(vidx, 1, Int(obs.nFFT))
  mapping =
    Utils.rescale_linear(beam.a_params[:peaks_idx], -vmax, vmax, 1, obs.nFFT)
  heights = obs.heights ./ 1000
  lo, hi = minmax(vidx...)
  raw = beam.spectra[:, lo:hi]
  den = beam.denoised_spectra[:, lo:hi]
  computed = mapping
  reference = beam.ref_M1 .* obs.DBS_Factor
  stw = beam.a_params[:mccf_ST_W]
  p1 = Plt.plot_doppler_spectra(
    spectra = 20 .* log10.(abs.(raw)),
    heights = heights,
    x_max = vl[2],
  )
  p2 = Plt.plot_doppler_spectra(
    spectra = 20 .* log10.(abs.(den)),
    comp_m = -computed,
    ref_m = -reference,
    heights = heights,
    x_max = vl[2],
  )
  p3 = Plt.plot_doppler_spectra(
    spectra = 20 .* log10.(abs.(den[10:60, :])),
    comp_m = -computed[10:60],
    ref_m = -reference[10:60],
    heights = heights[10:60],
    x_max = vl[2],
  )
  Plt.plot_ST_window(stw, vmax, heights, "white"; plot_handle = p3)
  display(plot(p1, p2, p3; layout = (1, 3), size = (1500, 600)))
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
  run_plot_STwindow()
end

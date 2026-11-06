using Radar
using Plots

function run_ime()
  filepath = joinpath("Data", "other", "EXP_DBS_CH4_29Jul2026_19_17_20")
  obs = get_observation(filepath)
  apply_beam_fx(obs, IME.compute_beam_spectra)
  apply_beam_fx(obs, IME.denoise_beam)
  apply_beam_fx(obs, IME.estimate_params; vel_axis = obs.velocities)
  Utils.compute_uvw(obs)

  panels = [
    Plt.plot_doppler_spectra(
      spectra   = 10 .* log10.(abs.(getproperty(obs, Symbol(d)).spectra) .^ 2),
      direction = d,
      comp_m    = -getproperty(obs, Symbol(d)).M1 .* obs.DBS_Factor,
      ref_m     = -getproperty(obs, Symbol(d)).ref_M1 .* obs.DBS_Factor,
      heights   = obs.heights,
      x_max     = obs.v_max,
    ) for d in Data.directions_str
  ]
  display(plot(panels...; layout = (1, 5), size = (1500, 500)))

  p = plot(layout = (1, 3), size = (1100, 700))
  Plt.plot_compared_vel(
    plotidx     = 1,
    calc        = obs.U,
    ref         = obs.ref_U,
    heights     = obs.heights,
    title       = "Zonal (U)",
    plot_handle = p,
  )
  Plt.plot_compared_vel(
    calc        = obs.V,
    ref         = obs.ref_V,
    plotidx     = 2,
    heights     = obs.heights,
    title       = "Meridional (V)",
    plot_handle = p,
  )
  Plt.plot_compared_vel(
    calc        = obs.W,
    ref         = obs.ref_W,
    plotidx     = 3,
    heights     = obs.heights,
    title       = "Vertical (W)",
    plot_handle = p,
  )
  display(p)
  println("IME processing finished")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
  run_ime()
end

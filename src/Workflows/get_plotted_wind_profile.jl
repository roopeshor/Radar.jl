function get_plotted_wind_profile(;
  nRangeBins::Integer = 173,
  start_height::Real = 0.0,
  end_height::Real = 20.0,
  make_plot::Bool = true,
)
  h = collect(
    range(Float64(start_height), Float64(end_height), length = nRangeBins),
  )
  ubg = Synth.synoptic_bg(h)
  ullj, vllj = Synth.somali_jet(h)
  utej, vtej = Synth.tropical_easterly_jet(h)
  uqbo, vqbo = Synth.stratospheric_qbo(h)
  umean = ubg .+ ullj .+ utej .+ uqbo
  vmean = vllj .+ vtej .+ vqbo
  woro = Synth.orographic_updraft(h, umean)
  ugw, vgw, wgw = Synth.gravity_waves(h)
  wplumes = Synth.convective_plumes(h)
  uturb, vturb, wturb = Synth.kolmogorov_turbulence(h, umean, vmean)
  u = umean .+ ugw .+ uturb
  v = vmean .+ vgw .+ vturb
  w = woro .+ wgw .+ wplumes .+ wturb

  if make_plot
    components = [
      (ubg, zeros(nRangeBins), "Synoptic background"),
      (ullj, vllj, "Low-level Somali jet"),
      (utej, vtej, "Tropical easterly jet"),
      (uqbo, vqbo, "Stratospheric QBO"),
      (woro, zeros(nRangeBins), "Orographic updraft"),
      (wplumes, zeros(nRangeBins), "Convective plumes"),
      (ugw, vgw, "Gravity waves"),
      (uturb, vturb, "Kolmogorov turbulence"),
      (u, v, "Total horizontal wind"),
    ]
    plots = [
      plot(
        c[1],
        h;
        label = "component 1",
        title = c[3],
        xlabel = "velocity (m/s)",
        ylabel = "height (km)",
      ) for c in components
    ]
    for (p, c) in zip(plots, components)
      plot!(p, c[2], h; label = "component 2")
    end
    display(plot(plots...; layout = (3, 3), size = (900, 1000)))
  end
  return (u, v, w)
end

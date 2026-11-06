function generate_wind_profile(;
  nRangeBins::Integer = 173,
  start_height::Real = 6.030,
  end_height::Real = 32.0,
  bg_z0_m::Real = 0.05,
  bg_ustar::Real = 0.3,
  bg_kappa::Real = 0.4,
  bg_blh::Real = 1.5,
  sj_w::Real = 12.0,
  sj_uv::Real = 14.0,
  sj_center::Real = 1.5,
  sj_s::Real = 0.8,
  tej_u::Real = 20.0,
  tej_v::Real = 3.0,
  tej_alt::Real = 15.0,
  tej_vthick::Real = 2.5,
  qbo_lambda::Real = 20.0,
  qbo_u::Real = 8.0,
  qbo_v::Real = 3.0,
  qbo_alt::Real = 17.0,
  ou_alpha_slope::Real = atan(3/20),
  ou_H_lift::Real = 2.0,
  gw_nmodes::Integer = 20,
  gw_amp::Real = 2.5,
  gw_max_lambda::Real = 4.0,
  gw_damp_scale::Real = 12.0,
  gw_cf::Real = 0.05,
  lw_Hs::Real = 7.0,
  jitter_amp::Real = 0.4,
  plume_count::Integer = 4,
  plume_max_alt::Real = 6.0,
  kt_vaf::Real = 0.25,
  kt_bf::Real = 1.0,
  kt_maxT::Real = 2.5,
  kt_Ric::Real = 0.25,
  kt_tlpeak::Real = 3.0,
  kt_vthick::Real = 2.5,
  lw_N_bv::Real = 0.012,
  turb_interm_sigma::Real = 0.5,
)
  h = collect(
    range(Float64(start_height), Float64(end_height), length = nRangeBins),
  )
  ubg = synoptic_bg(h; bg_z0_m, bg_ustar, bg_kappa, bg_blh)
  ullj, vllj = somali_jet(h; sj_w, sj_uv, sj_center, sj_s)
  utej, vtej = tropical_easterly_jet(h; tej_u, tej_v, tej_alt, tej_vthick)
  uqbo, vqbo = stratospheric_qbo(h; qbo_lambda, qbo_u, qbo_v, qbo_alt)
  umean = ubg .+ ullj .+ utej .+ uqbo
  vmean = vllj .+ vtej .+ vqbo
  umean .+= jitter_amp .* gaussian_smooth(randn(nRangeBins), 30)
  vmean .+= jitter_amp .* gaussian_smooth(randn(nRangeBins), 30)
  woro = orographic_updraft(h, umean; ou_alpha_slope, ou_H_lift)
  ugw, vgw, wgw = gravity_waves(
    h;
    gw_nmodes,
    gw_amp,
    gw_max_lambda,
    gw_damp_scale,
    gw_cf,
    lw_Hs,
  )
  wplumes = convective_plumes(h; plume_count, plume_max_alt)
  uturb, vturb, wturb = kolmogorov_turbulence(
    h,
    umean,
    vmean;
    kt_vaf,
    kt_bf,
    kt_maxT,
    kt_Ric,
    kt_tlpeak,
    kt_vthick,
    turb_interm_sigma,
    lw_N_bv,
  )
  u = umean .+ ugw .+ uturb
  v = vmean .+ vgw .+ vturb
  w = woro .+ wgw .+ wplumes .+ wturb
  return (h, u, v, w)
end

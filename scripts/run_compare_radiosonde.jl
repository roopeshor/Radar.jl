using Radar
using Plots
using Statistics

function compare_radiosonde(
  sonde_file = "./Data/Data_SOE/20_Jan_2026/Pdt_20Jan2026.txt",
  radar_raw_file = "./Data/Data_SOE/Radar_data/20_Jan_2026/EXP_DBS_CH4_20Jan2026_14_22_11.raw",
)
  obs_base = get_observation(radar_raw_file)
  sonde = RadarIO.read_radiosonde(sonde_file)
  obs_mccf = deepcopy(obs_base)
  apply_beam_fx(obs_mccf, MCCF.compute_all_params)
  obs_ime = deepcopy(obs_base)
  apply_beam_fx(obs_ime, IME.compute_beam_spectra)
  apply_beam_fx(obs_ime, IME.denoise_beam)
  apply_beam_fx(obs_ime, IME.estimate_params; vel_axis = obs_base.velocities)
  vertical_speed = [0.0; diff(sonde.altitude_m) ./ diff(sonde.time)]
  vertical_speed .-= mean(vertical_speed)
  radial = Utils.compute_radial_velocity(
    sonde.u,
    sonde.v,
    zeros(length(sonde.u)),
    obs_ime.north.m_fAzimuth,
    obs_ime.north.m_fOffZenith,
  )
  p1 = Plt.plot_doppler_spectra(
    spectra = log.(abs.(obs_mccf.north.denoised_spectra)),
    heights = obs_mccf.heights,
    comp_m = -obs_mccf.north.M1 .* obs_mccf.DBS_Factor,
    ref_m = -obs_mccf.north.ref_M1 .* obs_mccf.DBS_Factor,
    x_max = obs_mccf.v_max,
  )
  plot!(p1, radial, sonde.altitude_m; color = :white, label = "radiosonde")
  p2 = Plt.plot_doppler_spectra(
    spectra = log.(abs.(obs_ime.north.spectra)),
    heights = obs_ime.heights,
    comp_m = -obs_ime.north.M1 .* obs_ime.DBS_Factor,
    ref_m = -obs_ime.north.ref_M1 .* obs_ime.DBS_Factor,
    x_max = obs_ime.v_max,
  )
  plot!(p2, radial, sonde.altitude_m; color = :white, label = "radiosonde")
  display(plot(p1, p2; layout = (1, 2)))
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
  compare_radiosonde()
end

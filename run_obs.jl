include("src/Radar.jl")
using .Radar: Observation, apply_beam_fx, compute_uvw, compute_all_params, estimate_params
using .plt: plot_doppler_spectra, plot_compared_vel
filepath = joinpath("Data", "other", "EXP_DBS_CH4_29Jul2026_19_22_59")
println("Processing file: ", filepath)


obs = Observation(filepath)
println("completed proc ")
obs = apply_beam_fx(obs, compute_all_params)
println("applied fx ")
# obs = apply_beam_fx(obs, estimate_params, vel_axis=obs.velocities)
obs = compute_uvw(obs)
println("applied uvw ")

N = obs.north
S = obs.south
E = obs.east
W = obs.west
V = obs.vertical

# Spectra Figure
dirs = [N, S, E, W, V]
spectra_plots = []
for d in dirs
	println("plotting ", d.direction)
    p = plot_doppler_spectra(
        spectra   = 10 .* log10.(abs.(d.denoised_spectra) .^ 2),
        direction = d.direction,
        heights   = obs.heights,
        x_max     = obs.v_max
    )
    push!(spectra_plots, p)
end
plot(spectra_plots..., layout=(1, 5))

# UVW Figure
plot_compared_vel(
    calc    = @. -(E.M1 - W.M1) * obs.DBS_Factor_H,
    ref     = obs.ref_U,
    plotidx = 1,
    heights = obs.heights,
    title   = "Zonal (U)"
)

plot_compared_vel(
    calc    = @. -(N.M1 - S.M1) * obs.DBS_Factor_H,
    ref     = obs.ref_V,
    plotidx = 2,
    heights = obs.heights,
    title   = "Meridional (V)"
)

c_th = cosd(N.m_fOffZenith)
sum_M1 = E.M1 .+ W.M1 .+ N.M1 .+ S.M1
calc_W = @. -obs.DBS_Factor * (c_th * sum_M1 + V.M1) / (4 * c_th^2 + 1)

plot_compared_vel(
    calc    = calc_W,
    ref     = obs.ref_W,
    plotidx = 3,
    heights = obs.heights,
    title   = "Vertical (W)"
)

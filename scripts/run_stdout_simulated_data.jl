using Radar
using Statistics
using Random
using Plots

const _WARMED_METHODS = Set{Symbol}()
const FILE_PATHS = (
    (:params, "params_file.txt"),
    (:mse, "mse.tsv"),
    (:corr, "correlation.tsv"),
    (:bench, "bench.tsv")
)

function metrics(reference, estimate)
    mask = isfinite.(reference) .& isfinite.(estimate)
    count(mask) >= 2 || return (mse=NaN, corr=NaN)
    a = reference[mask]
    b = estimate[mask]
    return (mse=mean((a .- b) .^ 2), corr=Utils.safe_corrcoef(a, b))
end

function make_random_observation()
    altitude_modes = ("lower", "middle", "upper")
    stripe_counts = Dict("lower" => 4, "middle" => 3, "upper" => 2)
    altitude_range = rand(altitude_modes)
    stripe_base = stripe_counts[altitude_range]

    # Randomized atmospheric model parameters from the MATLAB generator.
    params = (
        kt_tlpeak=rand() * 5 + 10,
        kt_vaf=rand() / 5 + 0.1,
        kt_bf=rand() * 5 + 1,
        kt_maxT=rand() * 2 + 0.5,
        sj_center=rand() * 2 + 1,
        sj_uv=(rand() - 0.5) * 10,
        tej_u=(rand() - 0.5) * 30,
        tej_v=(rand() - 0.5) * 30,
        tej_alt=rand() * 10 + 2,
        qbo_lambda=rand() * 10 + 5,
        qbo_u=(rand() - 0.5) * 20,
        qbo_v=(rand() - 0.5) * 20,
        gw_nmodes=round(Int, rand() * 20) + 5,
        gw_amp=rand() * 2 + 1,
        plume_count=round(Int, rand() * 5) + 1,
        plume_max_alt=rand() * 3 + 5,
        bg_z0_m=(rand() + 0.5) / 100,
        gw_damp_scale=rand() * 10 + 1,
        qbo_alt=rand() * 10 + 5,
        bg_blh=rand() + 1,
    )

    headers = [Synth.create_synthetic_beamheader(
        headerFields=(m_sCurrentBeamCnt=i-1,), altitudeRange=altitude_range) for i in 1:5]
    hstart = headers[1].m_fWindow1StartHeight / 1000
    hend = headers[1].m_fWindow1EndHeight / 1000
    n = Int(headers[1].nRangeBins)

    # Build multiple bright-band rows with independently randomized widths.
    stripe_count = max(stripe_base, round(Int, randn() * stripe_base + 1))
    stripe_gap = max(1, round(Int, (n - max(10, stripe_count * 2)) / stripe_count))
    rows = Vector{Vector{Float64}}()
    stripes_per_bin = max.(stripe_base - 1, round.(Int, abs.(rand(stripe_count)) .* 3))
    for i in 1:stripe_count
        start_idx = (i - 1) * stripe_gap + round(Int,
            min(stripe_gap - stripes_per_bin[i], rand() * stripe_gap / 3))
        widths = max.(5.0, abs.(randn(stripes_per_bin[i] + 1) .* 6 .+ 4))
        indices = start_idx .+ collect(1:(stripes_per_bin[i]+1))
        for (idx, width) in zip(indices, widths)
            push!(rows, [Float64(idx), width])
        end
    end
    bands = isempty(rows) ? zeros(Float64, 0, 2) : reduce(vcat, permutedims.(rows))

    wind_kwargs = merge((nRangeBins=n, start_height=hstart, end_height=hend), params)
    _, u, v, w = Synth.generate_wind_profile(; wind_kwargs...)
    gc_rel_power = rand() / 200 + 0.0005
    params = merge(params, (gc_rel_power=gc_rel_power, bands=bands,
        altitude_range=altitude_range))
    bright_bands = Dict(d => bands for d in ("north", "south", "east", "west"))
    obs = Synth.create_synthetic_observation(u, v, w;
        headerFields=headers,
        SNR=10 .^ (collect(range(50, -35, length=n)) ./ 10),
        dc_rel_power=0.2, gc_rel_power=gc_rel_power,
        bright_bands=bright_bands)
    return obs, params
end

function _process_method!(obs, method::Symbol)
    if method == :ime
        apply_beam_fx(obs, IME.compute_beam_spectra)
        apply_beam_fx(obs, IME.denoise_beam)
        apply_beam_fx(obs, IME.estimate_params; vel_axis=obs.velocities)
    elseif method == :mccf
        apply_beam_fx(obs, MCCF.compute_all_params; wind_shear=1, win_coeff=2)
    elseif method == :simple
        apply_beam_fx(obs, Simple.compute_beam_spectra)
        apply_beam_fx(obs, Simple.denoise_beam)
        apply_beam_fx(obs, Simple.compute_moments)
    elseif method == :mc_den_ime
        apply_beam_fx(obs, MCCF.compute_beam_spectra)
        apply_beam_fx(obs, IME.denoise_beam)
        apply_beam_fx(obs, IME.estimate_params; vel_axis=obs.velocities)
    elseif method == :mc_ime
        apply_beam_fx(obs, MCCF.compute_beam_spectra)
        apply_beam_fx(obs, IME.estimate_params; vel_axis=obs.velocities)
    else
        throw(ArgumentError("Unknown processing method: $method"))
    end
    return obs
end


function run_method(base, method::Symbol)
    # Exclude Julia's first-call compilation cost from the algorithm timings.
    if !(method in _WARMED_METHODS)
        _process_method!(deepcopy(base), method)
        push!(_WARMED_METHODS, method)
    end
    obs = deepcopy(base)
    elapsed = @elapsed _process_method!(obs, method)
    return obs, elapsed
end

function run_stdout_simulated_data(
    folder_name="outputs/exp0";
    exp_arr=1:2,
    save_data=true,
    save_plots=true
)
    print("Started\n")
    mkpath(folder_name)
    plots_folder = joinpath(folder_name, "plots")
    mkpath(plots_folder)
    paths = Dict(k=>joinpath(folder_name, v) for (k, v) in FILE_PATHS)
    methods = (:ime, :mccf, :simple, :mc_den_ime, :mc_ime)
    pf = open(paths[:params], "w")
    mf = open(paths[:mse], "w")
    cf = open(paths[:corr], "w")
    bf = open(paths[:bench], "w")
    if save_data
        println(mf, "experiment_id\tmethod\tdirection\tMSE")
        println(cf, "experiment_id\tmethod\tdirection\tCorrelation")
        println(bf, "experiment_id\tmethod\ttime_ms")
    end
    for exp_id in exp_arr
        println("== Experiment: $exp_id")
        base, params = make_random_observation()
        println(pf, "$exp_id: $params")
        for method in methods
            result, elapsed = run_method(base, method)
            save_data && println(bf, "$exp_id\t$method\t$(1000 * elapsed)")
            for direction in Data.directions
                beam = getproperty(result, direction)
                ref = beam.ref_M1
                est = beam.M1
                if isempty(est)
					print("No useful data computed")
                    continue
                end
                m = metrics(ref, est)
                if save_data
                    println(mf, "$exp_id\t$method\t$direction\t$(m.mse)")
                    println(cf, "$exp_id\t$method\t$direction\t$(m.corr)")
                end
                if save_plots
                    p = Plt.plot_doppler_spectra(
                        spectra   = 20 .* log10.(max.(abs.(beam.spectra), eps())),
                        direction = beam.direction,
                        ref_m     = (-beam.ref_M1 .* result.DBS_Factor),
                        comp_m    = (-beam.M1 .* result.DBS_Factor),
                        heights   = result.heights ./ 1000,
                        x_max     = result.v_max
                    )
                    savefig(p, joinpath(plots_folder, "$(exp_id)_$(method)_$(direction).png"))
                end
            end
            println("$method: $(round(1000 * elapsed,digits=2)) ms")
        end
    end
    return nothing
end
if abspath(PROGRAM_FILE) == @__FILE__
    run_stdout_simulated_data("outputs/tester"; exp_arr=1:500)
end

"""
    stacked_spectrogram(spectrum; y_ticks, x_ticks, amplitude_scale=1.2, plot_color="blue", plot_gradiated=false, plot_add_colorbar=false, plot_thick=1, baseline_color=[0.5, 0.5, 0.5], baseline_thick=0, peak_marker="blue*", marker_size=4, ref_line_pos=0, ref_line_thick=1, ref_line_color="red", decimate_factor=1, x_axis_type="velocity")

plt.stacked_spectrogram  creates a stacked spectrogram.
    It can render plots with single color or gradient according to function value
    Optionally it can also mark max peak in each graph, and a reference line at given point

# Arguments
- `spectrum::AbstractMatrix{<:Real}`
- `y_ticks::AbstractVector{<:Real}`: y ticks vector
- `x_ticks::AbstractVector{<:Real}`: x ticks vector
- `amplitude_scale::Real=1.2`: scaling of amplitude of the spectral peaks.
- `plot_color::Any="blue"`: color of stacked plots
- `plot_gradiated::Bool=false`: whether to plot the lines with some gradient. if so the particular gradient should be specified in plot_color
- `plot_add_colorbar::Bool=false`: whether to add colorbar when plotting with amplitude gradient
- `plot_thick::Real=1`: thickness of stacked plot
- `baseline_color::Any=[0.5, 0.5, 0.5]`: color of baseline
- `baseline_thick::Real=0`: set this to 0 to remove line below every stacked plot
- `peak_marker::Any="blue*"`: peak marker, setting this to empty string will not draw any marker
- `marker_size::Real=4`: size of peak marker, setting this to 0 will not draw any marker
- `ref_line_pos::Real=0`: Where will be the reference line. By default its x=0
- `ref_line_thick::Real=1`: thickness of x=0 reference line. set this to 0 to remove reference line
- `ref_line_color::Any="red"`: color of reference line
- `decimate_factor::Integer=1`: fraction of input data to skip, 1 => nothing is skipped, increasing this decreases number of points plotted improving performance
- `x_axis_type::AbstractString="velocity"`: X-axis domain mode ("velocity" or "freq"). if its velocity, given input will be flipped
"""
function stacked_spectrogram(spectrum::AbstractMatrix{<:Real};
    y_ticks::AbstractVector{<:Real}=range(0, 1, length=size(spectrum, 1)) |> collect,
    x_ticks::AbstractVector{<:Real}=range(0, 1, length=size(spectrum, 2)) |> collect,
    amplitude_scale::Real=1.2,
    plot_color::Any=:blue,
    plot_gradiated::Bool=false,
    plot_add_colorbar::Bool=false,
    plot_thick::Real=1,
    baseline_color::Any=[0.5, 0.5, 0.5],
    baseline_thick::Real=0,
    peak_marker::Any="blue*",
    marker_size::Real=4,
    ref_line_pos::Real=0,
    ref_line_thick::Real=1,
    ref_line_color::Any=:red,
    decimate_factor::Integer=1,
    x_axis_type::AbstractString="velocity"
)

    if decimate_factor > 1
        # Decimate the horizontal data before constructing the blocks
        x_ticks = x_ticks[1:decimate_factor:end]
        spectrum = spectrum[:, 1:decimate_factor:end]
    end

    if x_axis_type == "velocity"
        spectrum = reverse(spectrum, dims=2)
    end

    # Amplitude scaling.
    # This finds the correct scale so that a waveform entirely fits in its given y range. amplitude_scaling does further scaling
    plot_dy = y_ticks[2] - y_ticks[1] # assume uniform tick size
    spectrum[isinf.(spectrum)] .= NaN

    y_count = size(spectrum, 1)
    y_max = Vector{Float64}(undef, y_count)
    y_max_idx = Vector{Int}(undef, y_count)
    y_min = Vector{Float64}(undef, y_count)

    for i in 1:y_count
        row = spectrum[i, :]
        valid = .!isnan.(row)
        if any(valid)
            ym, yi = findmax(row[valid])
            indices = findall(valid)
            y_max[i] = ym
            y_max_idx[i] = indices[yi]
            y_min[i] = minimum(row[valid])
        else
            y_max[i] = NaN
            y_max_idx[i] = 1
            y_min[i] = NaN
        end
    end

    y_diff_max = maximum(y_max .- y_min)
    amplitude_scale = amplitude_scale * plot_dy / y_diff_max

    y_stacked_all = y_ticks .+ (spectrum .* amplitude_scale)

    if baseline_thick > 0
        X_base = repeat([x_ticks[1], x_ticks[end], NaN], y_count)
        Y_base = vec(hcat(y_ticks', y_ticks', fill(NaN, 1, y_count))')
        plot!(X_base, Y_base, color=baseline_color, linewidth=baseline_thick, label=false)
    end

    # Vertical Reference Line at x = 0
    if ref_line_thick > 0
        vline!([ref_line_pos], color=ref_line_color, linewidth=ref_line_thick, label=false)
    end

    nans = fill(NaN, y_count, 1)
    # Create blocks with an extra NaN column at the end of each line
    # This forces plotting backends to drop the pen before starting the next line
    # also Flatten the blocks row-by-row into a single long row vector
    x_mat = repeat(x_ticks', y_count, 1)
    X_ticks_flat = vec(hcat(x_mat, nans)')
    Y_stacked_flat = vec(hcat(y_stacked_all, nans)')

    if plot_gradiated
        Y_vals = vec(hcat(spectrum, nans)')
        plot!(X_ticks_flat, Y_stacked_flat, line_z=Y_vals, color=plot_color, linewidth=plot_thick, label=false)
        if plot_add_colorbar
            plot!(colorbar=true)
        end
    else
        plot!(X_ticks_flat, Y_stacked_flat, color=plot_color, linewidth=plot_thick, label=false)
    end

    # Find and mark the peak
    if peak_marker != "" && marker_size > 0
        y_max_x = x_ticks[y_max_idx]
        y_max_offset = y_ticks .+ (y_max .* amplitude_scale)
        scatter!(y_max_x, y_max_offset, markersize=marker_size, label=false)
    end

    xlabel!("Doppler " * x_axis_type * " (m/s)")
    ylabel!("Altitude (km)")

    xlims!(minimum(x_ticks), maximum(x_ticks))
    ylims!(minimum(y_ticks) - amplitude_scale, maximum(y_ticks) + amplitude_scale)
end

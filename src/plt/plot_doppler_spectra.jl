"""
    plot_doppler_spectra(; spectra, comp_m=[], ref_m=[], direction="", heights=[], x_max=30, add_legend=true, comp_m_color=:black, ref_m_color=:red, plot_thickness=1, x_axis_type="velocity")

plt.plot_doppler_spectra Renders range-Doppler spectrographs with overlaid computed and reference moments.

# Arguments
- `spectra::AbstractMatrix{<:Real}`: Power spectra matrix (RangeBins x NFFT).
- `comp_m::AbstractVector{<:Real}=Float64[]`: Computed Doppler moments profile vector.
- `ref_m::AbstractVector{<:Real}=Float64[]`: Reference Doppler moments profile vector.
- `direction::AbstractString=""`: Beam direction label string.
- `heights::AbstractVector{<:Real}=Float64[]`: Altitude heights profile vector.
- `x_max::Real=30`: Maximum unambiguous Doppler velocity or frequency bound.
- `add_legend::Bool=true`: Flag to display plot legend.
- `comp_m_color::Any=:black`: Line color for computed moments curve.
- `ref_m_color::Any=:red`: Line color for reference moments curve.
- `plot_thickness::Real=1`: LineWidth of plots
- `x_axis_type::AbstractString="velocity"`: X-axis domain mode ("velocity" or "freq").
"""
function plot_doppler_spectra(; spectra::AbstractMatrix{<:Real},
    comp_m::AbstractVector{<:Real}=Float64[],
    ref_m::AbstractVector{<:Real}=Float64[],
    direction::AbstractString="",
    heights::AbstractVector{<:Real}=Float64[],
    x_max::Real=30,
    add_legend::Bool=true,
    comp_m_color::Any=:black,
    ref_m_color::Any=:red,
    plot_thickness::Real=1,
    x_axis_type::AbstractString="velocity"
)

    x_vals = range(-x_max, x_max, length=size(spectra, 2))
    y_vals = heights

    processed_spectra = spectra
    if x_axis_type == "velocity"
        processed_spectra = reverse(processed_spectra, dims=2)
    end

    heatmap(x_vals, y_vals, processed_spectra', color=:parula)

    if !isempty(ref_m)
        plot!(ref_m, heights, label="ref", linewidth=plot_thickness, color=ref_m_color)
    end
    if !isempty(comp_m)
        plot!(comp_m, heights, label="computed", linewidth=plot_thickness, color=comp_m_color)
    end

    xlims!(-x_max, x_max)
    ylims!(minimum(heights), maximum(heights))

    xlabel!("Doppler " * x_axis_type * " (m/s)")
    ylabel!("Altitude (km)")
    title!(direction)

    if !add_legend || (isempty(ref_m) && isempty(comp_m))
        plot!(legend=false)
    end
end

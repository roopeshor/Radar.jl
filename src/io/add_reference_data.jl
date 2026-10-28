using CSV
using DataFrames

# io.add_reference_data Attaches benchmark reference moments (.mmts) and UVW wind profiles (.uvw) to radar beams.
#
#   Parses reference text files matching the observation base filename, populating reference moment fields
#   (ref_M0, ref_M1, ref_M2, ref_SNR, ref_noise_level) on each beam and extracting reference wind vectors (U, V, W).
function add_reference_data(beams::Vector{RadarData}, filepath::String;
                            add_mmts::Bool = true,
                            add_uvw::Bool = true,
                            uvw_file_ext::String = "_W1.uvw",
                            mmts_file_ext = Data.dir2mmts_file_ext)

    beam_dict = Dict{String, RadarData}(
        "vertical" => RadarData(),
        "north"    => RadarData(),
        "west"     => RadarData(),
        "south"    => RadarData(),
        "east"     => RadarData()
    )

    folder, file_with_ext = splitdir(filepath)
    basename, _ = splitext(file_with_ext)

    height = Float64[]
    U = Float64[]
    V = Float64[]
    W = Float64[]

    for i in 1:length(beams)
        beam = beams[i]
        direction = beam.direction

        beam.ref_height = Float64[]
        beam.ref_M0 = Float64[]
        beam.ref_M1 = Float64[]
        beam.ref_M2 = Float64[]
        beam.ref_SNR = Float64[]
        beam.ref_noise_level = Float64[]

        if add_mmts
            ext = mmts_file_ext isa AbstractDict ? mmts_file_ext[direction] : getproperty(mmts_file_ext, Symbol(direction))
            path = joinpath(folder, basename * ext)
            if !isfile(path)
                @warn "$path not found, returning empty array"
            else
                output = CSV.read(path, DataFrame; delim='\t', header=1, silencewarnings=true)
                beam.ref_height      = Float64.(output[:, "Height (km)"])
                beam.ref_M0          = Float64.(output[:, "M0 Total Power(dBm)"])
                beam.ref_M1          = Float64.(output[:, "M1 Mean Doppler (Hz)"])
                beam.ref_M2          = Float64.(output[:, "M2 Doppler Spread (Hz)"])
                beam.ref_SNR         = Float64.(output[:, "SNR (dB)"])
                beam.ref_noise_level = Float64.(output[:, "Noise Level (dBm)"])
            end
        end
        beam_dict[direction] = beam
    end

    uvw_path = joinpath(folder, basename * uvw_file_ext)
    if add_uvw
        if !isfile(uvw_path)
            @warn "$uvw_path not found, returning empty array"
        else
            output = CSV.read(uvw_path, DataFrame; delim='\t', header=1, silencewarnings=true)
            height = Float64.(output[:, "Height (km)"])
            U      = Float64.(output[:, "Zonal(U) (m/s)"])
            V      = Float64.(output[:, "Meridional(V) (m/s)"])
            W      = Float64.(output[:, "Vertical (W) (m/s)"])
        end
    end

    return (
        beam_struct = (
            vertical = beam_dict["vertical"],
            north    = beam_dict["north"],
            west     = beam_dict["west"],
            south    = beam_dict["south"],
            east     = beam_dict["east"]
        ),
        U = U,
        V = V,
        W = W,
        height = height
    )
end

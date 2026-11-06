module RadarIO

using ..Types: RadarData, Data, BeamSet
using ..Utils: add_extra_RD_params!

export read_raw_file, write_raw_file, add_reference_data, read_radiosonde

include("binary_io.jl")
include("read_raw_file.jl")
include("write_raw_file.jl")
include("table_io.jl")
include("add_reference_data.jl")
include("read_radiosonde.jl")

end # module RadarIO

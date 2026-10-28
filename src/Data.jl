module Data

# Data : contains various mapping dictionaries. Putting everything here together so that
#   any changes made will be reflected everywhere
#
# cardinal direction, azimuth and beam number:
#
#           0°
#          N[4]
#
#
#  W[2]    V[0]     E[1]
#  270°              90°
#
#          S[3]
#          180°

# beam count to azimuth map
const cnt2az = Dict(
    0 => 0,
    1 => 270,
    2 => 90,
    3 => 180,
    4 => 0
)

# direction to beam count map
const dir2cnt = Dict(
    0 => "vertical",
    1 => "east",
    2 => "west",
    3 => "south",
    4 => "north"
)

# azimuth to direction map
const az2dir = Dict(
    0   => "north",
    90  => "east",
    180 => "south",
    270 => "west"
)

# direction to azimuth map
const dir2az = Dict(
    "north"    => 0,
    "vertical" => 0,
    "east"     => 90,
    "south"    => 180,
    "west"     => 270
)

# direction to mmts file extension map
const dir2mmts_file_ext = (
    vertical = "_Beam1_W1_Az_0.00_Oz_0.00.mmts",
    east     = "_Beam2_W1_Az_90.00_Oz_10.00.mmts",
    south    = "_Beam4_W1_Az_180.00_Oz_10.00.mmts",
    west     = "_Beam3_W1_Az_270.00_Oz_10.00.mmts",
    north    = "_Beam5_W1_Az_0.00_Oz_10.00.mmts"
)

# directions array
const directions = ["north", "east", "west", "south", "vertical"]

# azimuth and offzenith to direction map
function aoz2dir(; az::Real, oz::Real)
    if oz == 0
        return "vertical"
    end

    az_rounded = round(Int, az)
    return get(az2dir, az_rounded, "unknown")
end

end # module

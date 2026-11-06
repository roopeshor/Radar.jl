module Data

export cnt2az, cnt2dir, az2dir, dir2az, dir2mmts_file_ext, directions, aoz2dir

const cnt2az = Dict(0 => 0.0, 1 => 270.0, 2 => 90.0, 3 => 180.0, 4 => 0.0)
const cnt2dir =
  Dict(0 => "vertical", 1 => "west", 2 => "east", 3 => "south", 4 => "north")
const az2dir = Dict(0 => "north", 90 => "east", 180 => "south", 270 => "west")
const dir2az = Dict(
  "north"    => 0.0,
  "vertical" => 0.0,
  "east"     => 90.0,
  "south"    => 180.0,
  "west"     => 270.0,
)
const dir2mmts_file_ext = Dict(
  "vertical" => "_Beam1_W1_Az_0.00_Oz_0.00.mmts",
  "east"     => "_Beam2_W1_Az_90.00_Oz_10.00.mmts",
  "south"    => "_Beam4_W1_Az_180.00_Oz_10.00.mmts",
  "west"     => "_Beam3_W1_Az_270.00_Oz_10.00.mmts",
  "north"    => "_Beam5_W1_Az_0.00_Oz_10.00.mmts",
)
const directions = [
  Symbol("north"),
  Symbol("east"),
  Symbol("west"),
  Symbol("south"),
  Symbol("vertical"),
]
const directions_str = ["north", "east", "west", "south", "vertical"]

"Map an azimuth/off-zenith pair to a beam direction label."
function aoz2dir(; az::Real, oz::Real)
  iszero(oz) && return "vertical"
  return get(az2dir, round(Int, az), "unknown")
end

end

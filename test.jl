include("src/Radar.jl")
using .Radar: io, Observation
function main()
	x = Observation("Data/other/EXP_DBS_CH4_29Jul2026_19_22_59")
end

main();

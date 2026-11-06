# Radar.jl

This is a Julia port of original MATLAB implementation of [ST-radar processing](https://github.com/roopeshor/Radar-Signal-Processing). Got to this repo first to get more detailed readme
The port keeps processing functions in separate modules so that scripts are syntactically similar to MATLAB run scripts.

> This port just works, but needs to be polished, especially plotters.

## Layout

```text
.
├── Project.toml
├── Data/            # radar datas
├── src/
│   ├── Radar.jl     # package entry point
│   ├── Types/       # common datatyles
│   ├── RadarIO/     # reading .raw, .mmts, .uvw, radiosonde datas
│   ├── Utils/       # shared utilities
│   ├── IME/         # IME functions
│   ├── MCCF/        # MCCF functions
│   ├── Simple/      # simple processing functions
│   ├── Plt/         # various plotters
│   ├── Synth/       # synthetic data generation
│   └── Workflows/   # loading and doing stuff with observation
└── scripts/         # runnning scripts
```

## MATLAB-to-Julia mapping

The original `+io` package is named `RadarIO` because Julia already exports the `IO` module. The data classes are mutable structs because there is no classes in Julia as in MATLAB.

## Setup
Works in Julia 1.13. From the project directory run this once to install dependencies:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
## Run scripts

```sh
julia --project=. scripts/run_ime.jl
julia --project=. scripts/run_obs_data.jl
julia --project=. scripts/run_simulated_data.jl
```

For a raw observation named `EXP_NAME.raw`, the optional reference files are expected to be named as following:
```text
└── Data
    ├── EXP_NAME.raw
    ├── EXP_NAME_Beam1_W1_Az_0.00_Oz_0.00.mmts
    ├── EXP_NAME_Beam2_W1_Az_90.00_Oz_10.00.mmts
    ├── EXP_NAME_Beam3_W1_Az_270.00_Oz_10.00.mmts
    ├── EXP_NAME_Beam4_W1_Az_180.00_Oz_10.00.mmts
    ├── EXP_NAME_Beam5_W1_Az_0.00_Oz_10.00.mmts
    ├── EXP_NAME_W1.uvw

```

use `Revise.js` to use the package interactively

## Formatting
use [JuliaFormatter](https://juliaeditorsupport.github.io/JuliaFormatter.jl/).

install with package manager:
```
(Radar) pkg> app add JuliaFormatter
```
run with:
```jlfmt --inplace */**/*.jl```

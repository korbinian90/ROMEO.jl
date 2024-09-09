#!/usr/bin/env -S julia --color=yes --startup-file=no --threads=auto

## Usage

# Call with: `<path-to-file>/romeo.jl ARGS`
# On windows use: `julia --threads=auto <path-to-file>/romeo.jl ARGS`

# Example call:
# `./romeo.jl -p phase.nii.gz -m mag.nii.gz -t [0.025 0.05] -o output.nii.gz

import Pkg
Pkg.activate(@__DIR__)
try
    using ROMEO, MriResearchTools, ArgParse
catch
    Pkg.add(["ROMEO", "MriResearchTools", "ArgParse"])
    using ROMEO, MriResearchTools, ArgParse
end

@time msg = unwrapping_main(ARGS)
println(msg)

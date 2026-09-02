module ROMEO

using Statistics
using StatsBase

const NBINS = 256

# Evaluated while this package is precompiled, so the version is part of the
# image and does not depend on path metadata being readable at runtime. See
# `package_version` in provenance.jl.
const PKG_VERSION = pkgversion(@__MODULE__)

include("utility.jl")
include("priorityqueue.jl")
include("weights.jl")
include("seed.jl")
include("region_handling.jl")
include("algorithm.jl")
include("unwrapping.jl")
include("voxelquality.jl")
include("provenance.jl")
include("parse.jl")

unwrapping_main(args...; kwargs...) = @warn("Type `using ArgParse` to use this function \n `?unwrapping_main` for argument help")

export unwrap, unwrap!, unwrap_individual, unwrap_individual!, voxelquality, unwrapping_main,
       write_provenance, write_citations, register_citation!, package_version

end # module

module ROMEO

using Statistics
using StatsBase

const NBINS = 256

# Baked in at precompile time; include_dependency so a version bump invalidates the cache.
include_dependency(joinpath(@__DIR__, "..", "Project.toml"))
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

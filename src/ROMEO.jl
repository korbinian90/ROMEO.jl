module ROMEO

using Statistics
using StatsBase

const NBINS = 256

include("utility.jl")
include("priorityqueue.jl")
include("weights.jl")
include("seed.jl")
include("region_handling.jl")
include("algorithm.jl")
include("unwrapping.jl")
include("voxelquality.jl")

export unwrap, unwrap!, unwrap_individual, unwrap_individual!, voxelquality, unwrapping_main

unwrapping_main(args...; kwargs...) = @warn("Type `using MriResearchTools ArgParse` to use this function \n `?unwrapping_main` for argument help")

end # module

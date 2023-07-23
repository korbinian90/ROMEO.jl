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

unwrapping_main(args...; kwargs...) = @warn("Type `using ArgParse` to use this function \n `?unwrapping_main` for argument help")
if !isdefined(Base, :get_extension) # fallback for julia < 1.9
    include("../ext/RomeoApp/RomeoApp.jl")
end

export unwrap, unwrap!, unwrap_individual, unwrap_individual!, voxelquality, unwrapping_main

end # module

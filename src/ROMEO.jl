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
include("singularities.jl")
include("singularity_correction.jl")

unwrapping_main(args...; kwargs...) = @warn("Type `using ArgParse` to use this function \n `?unwrapping_main` for argument help")

export unwrap, unwrap!, unwrap_individual, unwrap_individual!, voxelquality, unwrapping_main
export detect_singularities, residues, singularity_mask, branchcuts, congruence_error
export fix_singularities, fix_singularities!

end # module

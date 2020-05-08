module ROMEO

using Statistics
using StatsBase

const NBINS = 256

include("utility.jl")
include("priorityqueue.jl")
include("dilation_and_erosion.jl")
include("weights.jl")
include("algorithm.jl")
include("seed.jl")
include("unwrapping.jl")

export unwrap, unwrap!, unwrap_individual, unwrap_individual!

end # module

module ROMEO

using Statistics
using StatsBase
using Images # for minfilter in wrapfit.jl

const NBINS = 256

include("utility.jl")
include("priorityqueue.jl")
include("weights.jl")
include("seed.jl")
include("region_handling.jl")
include("algorithm.jl")
include("unwrapping.jl")
include("wrapfit.jl")

export unwrap, unwrap!, unwrap_individual, unwrap_individual!

end # module

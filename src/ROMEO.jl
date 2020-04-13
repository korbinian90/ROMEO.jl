module ROMEO
using Statistics

struct Data
    mag::AbstractArray
    B0::AbstractArray
    PO::AbstractArray
    TEs::AbstractArray
end

include("utility.jl")
include("priorityqueue.jl")
include("dilation_and_erosion.jl")
include("weights.jl")
include("algorithm.jl")
include("seed.jl")
include("unwrapping.jl")

export unwrap, unwrap!, unwrap4d!

end # module

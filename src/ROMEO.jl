module ROMEO

using Statistics

const NBINS = 256

# Baked in at precompile time; include_dependency so a version bump invalidates
# the cache. Read from Project.toml rather than through pkgversion, which
# returns nothing when the package is precompiled with --strip-metadata, as
# juliac does.
const PKG_VERSION = let toml = joinpath(@__DIR__, "..", "Project.toml")
    include_dependency(toml)
    m = match(r"^version\s*=\s*\"([^\"]+)\""m, read(toml, String))
    m === nothing && error("no version field in $toml")
    VersionNumber(m.captures[1])
end

include("cli.jl")
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

unwrapping_main(args...; kwargs...) = @warn("Type `using MriResearchTools` to use this function \n `?unwrapping_main` for argument help")

export unwrap, unwrap!, unwrap_individual, unwrap_individual!, voxelquality, unwrapping_main,
       write_provenance, write_citations, register_citation!, package_version

end # module

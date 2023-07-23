using ROMEO
using Test
using MriResearchTools

nan_test(I1, I2) = I1[.!isnan.(I1)] ≈ I2[.!isnan.(I2)]

@testset "ROMEO.jl" begin
    include("features.jl")
    include("specialcases.jl")
    include("dsp_tests.jl")
    include("mri.jl")
    include("voxelquality.jl")
    #include("timing.jl")
end

@testset "RomeoApp" begin
    using ArgParse
    include("RomeoApp/dataset_small.jl")
    include("RomeoApp/dataset_small2.jl")
end

## print version to verify
println()
unwrapping_main(["--version"])

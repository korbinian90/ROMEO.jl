using ROMEO
using Test
using Aqua
using MriResearchTools

nan_test(I1, I2) = I1[.!isnan.(I1)] ≈ I2[.!isnan.(I2)]

@testset "ROMEO.jl" begin
    include("features.jl")
    include("specialcases.jl")
    include("dsp_tests.jl")
    include("mri.jl")
    include("voxelquality.jl")
    include("threading.jl")
    include("provenance.jl")
    include("parse.jl")
    include("cli.jl")
    #include("timing.jl")
end

if VERSION ≥ v"1.9"
    @testset "RomeoApp" begin
        include("RomeoApp/dataset_small.jl")
        include("RomeoApp/dataset_small2.jl")
    end
end

## print version to verify
println()
unwrapping_main(["--version"])

@testset "Aqua" begin
    Aqua.test_all(ROMEO)
end

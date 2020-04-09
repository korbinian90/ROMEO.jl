using ROMEO
using Test
using MriResearchTools

nan_test(I1, I2) = I1[.!isnan.(I1)] ≈ I2[.!isnan.(I2)]

@testset "ROMEO.jl" begin
    include("features.jl")
    include("specialcases.jl")
    include("dsp_tests.jl")
    include("mri.jl")
    include("timing.jl")
end

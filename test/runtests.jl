using ROMEO
using Test

nan_test(I1, I2) = @test I1[.!isnan.(I1)] ≈ I2[.!isnan.(I2)]

@testset "ROMEO.jl" begin
    include("specialcases.jl")
    include("dsp_tests.jl")
    include("mri.jl")
end

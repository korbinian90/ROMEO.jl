using ROMEO
using Test

@testset "ROMEO.jl" begin
    include("specialcases.jl")
    include("dsp_tests.jl")
    include("mri.jl")
    #include("timing.jl")
    # Write your own tests here.
end

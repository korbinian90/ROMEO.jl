using ROMEO
using Test

@testset "ROMEO.jl" begin
    include("unwrapscore_test.jl")
    include("unwrap_test.jl")
    include("unwrap_dsp_test.jl")
    #include("timing_test.jl")
    # Write your own tests here.
end

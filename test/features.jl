@testset "Features" begin

for l in 7:5:20
    for offset in 2π .* [0, 2, -1, 10, -50]
        phase_uw = collect(range(-2π; stop=2π, length=l))
        phase = rem2pi.(phase_uw, RoundNearest) .+ offset
        for mag in [nothing, collect(range(10; stop=1, length=l)), collect(range(1; stop=10, length=l))]
            unwrapped = unwrap(phase; mag=mag, correctglobal=true)
            @test unwrapped ≈ phase_uw
        end
    end
end

end

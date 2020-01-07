using MRI

phase = ones(3,3,3)
unwrap(phase)

point = ones(1,1,1)
@test_throws AssertionError unwrap(point)

phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phaseni = readphase(phasefile)
magni = readmag(magfile)

function unwrap_test(wrapped; keyargs...)
    unwrapped = unwrap(wrapped)
    # test that unwrapped is not a copy of phase
    @test unwrapped != wrapped
    # test that all resulting values are only 2π different
    @test all(isapprox.(rem2pi.(unwrapped - wrapped, RoundNearest), 0; atol=1e-6))
    unwrapped
end

phase = phaseni.raw

t1 = unwrap_test(phase)
t2 = unwrap_test(phase; mag=magni)
t3 = unwrap_test(phase; weights=:bestpath)

#@test t1 != t2
#@test t2 != t3
#@test t1 != t3

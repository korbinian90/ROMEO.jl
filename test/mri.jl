@testset "MRI tests" begin

using NIfTI

phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
echo = 3
phase = niread(phasefile).raw[:,:,:,echo]
mag = niread(magfile).raw[:,:,:,echo]

## accuracy tests
function unwrap_test(wrapped; keyargs...)
    unwrapped = unwrap(wrapped; keyargs...)
    # test that unwrapped is not a copy of phase
    @test unwrapped != wrapped
    # test that all resulting values are only 2π different
    @test all(isapprox.(rem2pi.(unwrapped - wrapped, RoundNearest), 0; atol=1e-6))
    unwrapped
end

t1 = unwrap_test(phase)
t2 = unwrap_test(phase; mag=mag)
t3 = unwrap_test(phase; weights=:bestpath)

# all results should be different
@test t1 != t2
@test t2 != t3
@test t1 != t3

## performance tests (at end to avoid first run overhead)
@test (@timed unwrap(phase))[5].poolalloc < 2e3
@test (@timed unwrap(phase; mag=mag))[5].poolalloc < 2e3
@test (@timed unwrap(phase; weights=:bestpath))[5].poolalloc < 3e3

## NaN tests
nanphase = copy(phase)
nanphase[1,:,:] .= NaN
nan_t1 = copy(t1)
nan_t1[1,:,:] .= NaN
nan_test(unwrap(nanphase), nan_t1)

nanmag = copy(mag)
nanmag[1,:,:] .= NaN
nan_test(unwrap(phase; mag=nanmag)[2:end,:,:], t2[2:end,:,:])

end

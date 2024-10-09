@testset "MRI tests" begin

phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
echo = 3
echo2 = 2
phase4D = niread(phasefile).raw
mag4D = niread(magfile).raw
phase = phase4D[:,:,:,echo]
mag = mag4D[:,:,:,echo]
phase2 = niread(phasefile).raw[:,:,:,echo2]
TEs = [echo, echo2]

## accuracy tests
function unwrap_test(wrapped; keyargs...)
    unwrapped = unwrap(wrapped; keyargs...)
    # test that unwrapped is not a copy of phase
    @test unwrapped != wrapped
    # test that all resulting values are only 2π different
    @test all(isapprox.(rem2pi.(unwrapped - wrapped, RoundNearest), 0; atol=1e-5))
    unwrapped
end

t = []
push!(t, unwrap_test(phase))
push!(t, unwrap_test(phase; mag=mag))
#push!(t, unwrap_test(phase4D)) #TODO use test data set with noise to see difference
push!(t, unwrap_test(phase4D; mag=mag4D, TEs=TEs))
push!(t, unwrap_individual(phase4D; mag=mag4D, TEs=TEs))
push!(t, unwrap_test(phase; weights=:bestpath))
#push!(t, unwrap_test(phase; weights=:romeo, mask=robustmask(mag)))
push!(t, unwrap_test(phase; weights=:romeo2, mag=mag, TEs=TEs, phase2=phase2))
push!(t, unwrap_test(phase; weights=:romeo3, mag=mag, TEs=TEs, phase2=phase2))
push!(t, unwrap_test(phase; weights=:romeo4, mag=mag, TEs=TEs, phase2=phase2))
push!(t, unwrap_test(phase; mag=mag, maxseeds=50))
push!(t, unwrap_test(phase4D; mag=mag4D, TEs=TEs, template=2))
push!(t, unwrap_test(phase4D; mag=mag4D, TEs=TEs, template=2, temporal_uncertain_unwrapping=0.9))

# all results should be different
for i in 1:length(t), j in 1:(i-1)
    if (t[i] == t[j])
        @show i j
    end
    @test t[i] != t[j]
end

t = []
push!(t, unwrap_test(phase; mag=mag, maxseeds=50))
push!(t, unwrap_test(phase; mag=mag, maxseeds=50, merge_regions=true, correct_regions=true))
#TODO correct_regions does not affect outcome
#push!(t, unwrap_test(phase; mag=mag, maxseeds=50, merge_regions=true))
for i in 1:length(t), j in 1:(i-1)
    @test t[i] != t[j]
end 

@test unwrap_individual(phase4D; mag=mag4D, TEs=TEs) == unwrap(phase4D; mag=mag4D, TEs=TEs, individual=true)

## performance tests (not at beginning to avoid first run overhead)
if VERSION ≥ v"1.8" # different performance on older julia versions
    @test (@timed unwrap(phase))[5].poolalloc < 6e3
    @test (@timed unwrap(phase; mag=mag))[5].poolalloc < 7e3
    @test (@timed unwrap(phase; weights=:bestpath))[5].poolalloc < 3e4
end

## NaN tests
nanphase = copy(phase)
nanphase[1,:,:] .= NaN
nan_unwrapped = unwrap_test(phase)
nan_unwrapped[1,:,:] .= NaN
@test nan_test(unwrap(nanphase), nan_unwrapped)

nanmag = copy(mag)
nanmag[1,:,:] .= NaN
@test nan_test(unwrap(phase; mag=nanmag)[2:end,:,:], unwrap_test(phase; mag=mag)[2:end,:,:])

end

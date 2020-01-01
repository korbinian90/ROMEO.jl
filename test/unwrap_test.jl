using MRI

phasefile = joinpath("data", "small", "Phase.nii")
phase = readphase(phasefile)

unwrapped = romeo(phase)

# test that unwrapped is not a copy of phase
@test unwrapped != phase

# test that all resulting values are only 2π different
@test all(isapprox.(rem2pi.(unwrapped - phase, RoundNearest), 0; atol=1e-6))

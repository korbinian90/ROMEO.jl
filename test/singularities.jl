using Statistics
include("vein_phantom.jl")

@testset "Singularities" begin

# S = (ρ - R) + iz has its zero line on the ring ρ=R, z=0, so the phase of S
# contains exactly one closed vortex line of a known size.
function vortex_ring(sz=(31, 31, 15); R=8.0)
    c = (sz .+ 1) ./ 2
    return [atan(I[3] - c[3], sqrt((I[1] - c[1])^2 + (I[2] - c[2])^2) - R) for I in CartesianIndices(sz)]
end

ring = vortex_ring()
sz = size(ring)
s = detect_singularities(ring)
ph = vein_phantom()

## detection
@test length(s.loops) == 1
@test s.loops[1].closed
@test ROMEO.loop_length(s.loops[1]) == 60
@test sum(abs, s.residues) == 60
@test count(singularity_mask(s)) == length(s.loops[1].voxels)
@test !any(singularity_mask(s; min_length=100))

## the residue field is divergence free, therefore the lines are closed
@test all(ROMEO.cube_divergence(s.residues, J) == 0 for J in CartesianIndices(ntuple(d -> 1:sz[d]-1, 3)))

## residues only depend on the wrapped differences (congruence invariance)
unwrapped = unwrap(ring)
@test residues(unwrapped) == residues(ring)
# the ring has phase differences of exactly π (the ambiguous case of rem2pi),
# the phantom phase is generic
congruent = ph.phase .+ 2π .* [iseven(sum(Tuple(I))) ? 1 : -2 for I in CartesianIndices(size(ph.phase))]
@test residues(congruent) == residues(ph.phase)
@test residues(unwrap(ph.phase; mag=ph.mag)) == residues(ph.phase)

## masked detection
mask = falses(sz)
mask[10:20, 10:20, :] .= true
@test sum(abs, residues(ring; mask)) < sum(abs, residues(ring))

## smooth phase: no residues, no branch cuts, correction is a no-op
ramp = [0.3I[1] + 0.2I[2] - 0.1I[3] for I in CartesianIndices((20, 20, 20))]
smooth_uw = unwrap(rem2pi.(ramp, RoundNearest))
@test sum(abs, residues(smooth_uw)) == 0
@test count(branchcuts(smooth_uw)) == 0
before = copy(smooth_uw)
info = fix_singularities!(smooth_uw; mode=:cascade)
@test smooth_uw == before # bit identical
@test info.nchanged == 0
@test !any(info.singularities)

## a residue leaves a branch cut in the congruent ROMEO result
@test count(branchcuts(unwrapped)) > 0

## correction
for mode in (:lsq, :inpaint, :smooth, :cascade)
    uw = unwrap(ring)
    orig = copy(uw)
    info = fix_singularities!(uw; mode)
    @test info.nfixed == 1
    @test info.nskipped == 0
    @test count(branchcuts(uw)) == 0 # the branch cut is gone
    @test uw[.!info.changed] == orig[.!info.changed] # bit identical outside
    @test all(info.delta[.!info.changed] .== 0)
    @test 0 < info.nchanged < length(uw) / 3
end

## :mask only detects, :off does nothing at all
uw = unwrap(ring)
orig = copy(uw)
info = fix_singularities!(uw; mode=:mask)
@test uw == orig
@test any(info.singularities)
@test any(info.cuts)
@test isempty(info.delta)
@test info.nchanged == 0

info = fix_singularities!(uw; mode=:off)
@test uw == orig
@test info.nloops == 0
@test !any(info.singularities)

@test_throws ArgumentError fix_singularities!(copy(uw); mode=:nonsense)

## fallback rules: too long, too large or not dark enough is only marked
info = fix_singularities!(copy(uw); mode=:lsq, max_loop_length=10)
@test info.nfixed == 0 && info.nskipped == 1
info = fix_singularities!(copy(uw); mode=:lsq, max_patch_voxels=10)
@test info.nfixed == 0 && info.nskipped == 1
info = fix_singularities!(copy(uw); mode=:lsq, mag=ones(sz)) # no signal void -> don't touch
@test info.nfixed == 0 && info.nskipped == 1
info = fix_singularities!(copy(uw); mode=:lsq, mag=ones(sz), min_void_voxels=0)
@test info.nfixed == 1

## non-mutating version
fixed, info = fix_singularities(uw; mode=:cascade)
@test uw == orig
@test fixed != uw
@test count(branchcuts(fixed)) == 0

## NaN input must not crash
nanring = copy(uw)
nanring[1, :, :] .= NaN
info = fix_singularities!(nanring; mode=:cascade)
@test isfinite(sum(filter(isfinite, nanring)))

## vein phantom: the singularities appear on their own, not put in by hand
TE = 25e-3
truth = 2π .* ph.field .* TE
@test sum(abs, residues(ph.phase)) > 0

uw = unwrap(ph.phase; mag=ph.mag)
orig = copy(uw)
info = fix_singularities!(uw; mag=ph.mag, mode=:cascade)
@test info.nfixed > 0
@test sum(abs, residues(uw)) < 0.1 * sum(abs, residues(orig))
@test count(branchcuts(uw)) < 0.1 * count(branchcuts(orig))

# protection metrics
@test uw[.!info.changed] == orig[.!info.changed]
@test info.nchanged == count(info.changed)
mag_median = median(ph.mag)
@test count(info.changed[ph.mag.≥mag_median]) == 0 # only low signal voxels are touched
top_decile = ph.mag .≥ quantile(vec(ph.mag), 0.9)
@test maximum(congruence_error(uw, orig)[top_decile]) == 0

# the correction has to move the phase towards the true field
align(x) = x .- 2π * round(median(x .- truth) / 2π)
rmse(x) = sqrt(sum(abs2, align(x) .- truth) / length(x))
@test rmse(uw) < 0.7 * rmse(orig)

## negative control: 1mm, 3T, short TE should not produce singularities at all
lowres = vein_phantom(; sz=(24, 24, 12), voxelsize=1.0, radius=0.5, B0=3.0, TE=8e-3)
@test sum(abs, residues(lowres.phase)) == 0

## 4D
ph4 = cat(ph.phase, ph.phase; dims=4)
mag4 = cat(ph.mag, ph.mag; dims=4)
uw4 = unwrap(ph4; mag=mag4, TEs=[1, 2])
orig4 = copy(uw4)
info4 = fix_singularities!(uw4; mag=mag4, mode=:mask)
@test size(info4.singularities) == size(ph4)
@test info4.nchanged == 0 && uw4 == orig4
info4 = fix_singularities!(uw4; mag=mag4, mode=:cascade)
@test size(info4.delta) == size(ph4)
@test uw4[.!info4.changed] == orig4[.!info4.changed]
@test info4.nfixed > 0

end

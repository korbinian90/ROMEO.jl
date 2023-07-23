@testset "Dataset small2" begin

p = joinpath("data", "small2")
phasefile_me = joinpath(p, "Phase.nii")
fn_mag = joinpath(p, "Mag.nii")
fn_phase = joinpath(p, "Phase.nii")
tmpdir = mktempdir()

unwrapping_main(["-p", fn_phase, "-m", fn_mag, "-o", tmpdir, "-v"])

end
@testset "ROMEO function tests" begin

original_path = abspath(".")
p = abspath(joinpath("data", "small"))
tmpdir = mktempdir()
cd(tmpdir)
phasefile_me = joinpath(p, "Phase.nii")
phasefile_me_nan = joinpath(p, "phase_with_nan.nii")
magfile_me = joinpath(p, "Mag.nii")
phasefile_1eco = joinpath(tmpdir, "Phase.nii")
phasefile_2D = joinpath(tmpdir, "Phase2D.nii")
magfile_1eco = joinpath(tmpdir, "Mag.nii")
magfile_2D = joinpath(tmpdir, "Mag2D.nii")
phasefile_1arreco = joinpath(tmpdir, "Phase.nii")
magfile_1arreco = joinpath(tmpdir, "Mag.nii")
maskfile = joinpath(tmpdir, "Mask.nii")
savenii(niread(magfile_me)[:,:,:,1] |> I -> I .> MriResearchTools.median(I), maskfile)
savenii(niread(phasefile_me)[:,:,:,1], phasefile_1eco)
savenii(niread(magfile_me)[:,:,:,1], magfile_1eco)
savenii(niread(phasefile_me)[:,:,:,[1]], phasefile_1arreco)
savenii(niread(magfile_me)[:,:,:,[1]], magfile_1arreco)
savenii(niread(phasefile_me)[:,:,1,1], phasefile_2D)
savenii(niread(magfile_me)[:,:,1,1], magfile_2D)

phasefile_me_5D = joinpath(tmpdir, "phase_multi_channel.nii")
magfile_5D = joinpath(tmpdir, "mag_multi_channel.nii")
savenii(repeat(niread(phasefile_me),1,1,1,1,2), phasefile_me_5D)
savenii(repeat(niread(magfile_me),1,1,1,1,2), magfile_5D)

function test_romeo(args)
    folder = tempname()
    args = [args..., "-o", folder]
    try
        msg = unwrapping_main(args)
        @test msg == 0
        @test isfile(joinpath(folder, "unwrapped.nii"))
    catch e
        println(args)
        println(sprint(showerror, e, catch_backtrace()))
        @test "test failed" == "with error" # signal a failed test
    end
end

configurations_se(pf, mf) = vcat(configurations_se.([[pf], [pf, "-m", mf]])...)
configurations_se(pm) = [
    [pm...],
    [pm..., "-g"],
    [pm..., "-N"],
    [pm..., "-i"],
    [pm..., "-q"],
    [pm..., "-Q"],
    [pm..., "-u"],
    [pm..., "-w", "romeo"],
    [pm..., "-w", "bestpath"],
    [pm..., "-w", "1010"],
    [pm..., "-w", "101011"],
    [pm..., "--threshold", "4"],
    [pm..., "-s", "50"],
    [pm..., "-s", "50", "--merge-regions"],
    [pm..., "-s", "50", "--merge-regions", "--correct-regions"],
    [pm..., "--wrap-addition", "0.1"],
    [pm..., "-k", "robustmask"],
    [pm..., "-k", "nomask"],
    [pm..., "-k", "qualitymask"],
    [pm..., "-k", "qualitymask", "0.1"],
]
configurations_me(phasefile_me, magfile_me) = vcat(configurations_me.([[phasefile_me], [phasefile_me, "-m", magfile_me]])...)
configurations_me(pm) = [
    [pm..., "-e", "1:2", "-t", "[2,4]"], # giving two echo times for two echoes used out of three
    [pm..., "-e", "[1,3]", "-t", "[2,4,6]"], # giving three echo times for two echoes used out of three
    [pm..., "-e", "[1", "3]", "-t", "[2,4,6]"],
    [pm..., "-t", "[2,4,6]"],
    [pm..., "-t", "2:2:6"],
    [pm..., "-t", "[2.1,4.2,6.3]"],
    [pm..., "-t", "epi"], # shorthand for "ones(<num-echoes>)"
    [pm..., "-t", "epi", "5.3"], # shorthand for "5.3*ones(<num-echoes>)"
    [pm..., "-B", "-t", "[2,4,6]"],
    [pm..., "-B", "-t", "[2" ,"4", "6]"], # when written like [2 4 6] in command line
    [pm..., "--temporal-uncertain-unwrapping", "-t", "[2,4,6]"],
    [pm..., "--template", "1", "-t", "[2,4,6]"],
    [pm..., "--template", "3", "-t", "[2,4,6]"],
    [pm..., "--phase-offset-correction", "-t", "[2,4,6]"],
    [pm..., "--phase-offset-correction", "bipolar", "-t", "[2,4,6]"],
    [pm..., "--phase-offset-correction", "-t", "[2,4,6]", "--phase-offset-smoothing-sigma-mm", "[5,8,4]"],
    [pm..., "--phase-offset-correction", "-t", "[2,4,6]", "--write-phase-offsets"],
]

files = [(phasefile_1eco, magfile_1eco), (phasefile_1arreco, magfile_1arreco), (phasefile_1eco, magfile_1arreco), (phasefile_1arreco, magfile_1eco), (phasefile_2D, magfile_2D)]
for (pf, mf) in files, args in configurations_se(pf, mf)
    test_romeo(args)
end
for args in configurations_me(phasefile_me, magfile_me)
    test_romeo(args)
end
for args in configurations_se(["-p", phasefile_me, "-m", magfile_me, "-t", "[2,4,6]"])
    test_romeo(args)
end
for args in configurations_me(phasefile_me_5D, magfile_5D)[end-2:end] # test the last 3 configurations_me lines for coil combination
    test_romeo(args)
end
files_se = [(phasefile_1eco, magfile_1eco), (phasefile_1arreco, magfile_1arreco)]
for (pf, mf) in files_se
    b_args = ["-B", "-t", "3.06"]
    test_romeo(["-p", pf, b_args...])
    test_romeo(["-p", pf, "-m", mf, b_args...])
end

test_romeo([phasefile_me_nan, "-t", "[2,4]", "-k", "nomask"])

## Test error and warning messages
m = "multi-echo data is used, but no echo times are given. Please specify the echo times using the -t option."
@test_throws ErrorException(m) unwrapping_main(["-p", phasefile_me, "-o", tmpdir, "-v"])

m = "masking option '0.8' is undefined (Maybe '-k qualitymask 0.8' was meant?)"
@test_throws ErrorException(m) unwrapping_main(["-p", phasefile_1eco, "-o", tmpdir, "-v", "-k", "0.8"])

m = "masking option 'blub' is undefined"
@test_throws ErrorException(m) unwrapping_main(["-p", phasefile_1eco, "-o", tmpdir, "-v", "-k", "blub"])

m = "Phase offset determination requires all echo times! (2 given, 3 required)"
@test_throws ErrorException(m) unwrapping_main(["-p", phasefile_me_5D, "-o", tmpdir, "-v", "-t", "[1,2]", "-e", "[1,2]", "--phase-offset-correction"])

m = "echoes=[1,5]: specified echo out of range! Number of echoes is 3"
@test_throws ErrorException(m) unwrapping_main(["-p", phasefile_me, "-o", tmpdir, "-v", "-t", "[1,2,3]", "-e", "[1,5]"])

m = "echoes=[1,5} wrongly formatted!"
@test_throws ErrorException(m) unwrapping_main(["-p", phasefile_me, "-o", tmpdir, "-v", "-t", "[1,2,3]", "-e", "[1,5}"])

m = "Number of chosen echoes is 2 (3 in .nii data), but 5 TEs were specified!"
@test_throws ErrorException(m) unwrapping_main(["-p", phasefile_me, "-o", tmpdir, "-v", "-t", "[1,2,3,4,5]", "-e", "[1,2]"])

m = "size of magnitude and phase does not match!"
@test_throws ErrorException(m) unwrapping_main(["-p", phasefile_me, "-o", tmpdir, "-v", "-t", "[1,2,3]", "-m", magfile_1eco])

m = "robustmask was chosen but no magnitude is available. No mask is used!"
@test_logs (:warn, m) match_mode=:any unwrapping_main(["-p", phasefile_1eco, "-o", tmpdir])

m = "The echo times 1 and 2 ([1.1, 1.1, 1.1]) need to be different for MCPC-3D-S phase offset correction! No phase offset correction performed"
@test_logs (:warn, m) match_mode=:any unwrapping_main(["-p", phasefile_me, "-m", magfile_me, "-o", tmpdir, "-t", "[1.1, 1.1, 1.1]", "--phase-offset-correction"])

@test_logs unwrapping_main(["-p", phasefile_1eco, "-o", tmpdir, "-m", magfile_1eco]) # test that no warning appears

## test maskfile
unwrapping_main([phasefile_1eco, "-k", maskfile])

## test no-rescale
phasefile_me_uw = joinpath(tempname(), "unwrapped.nii")
phasefile_me_uw_wrong = joinpath(tempname(), "wrong_unwrapped.nii")
phasefile_me_uw_again = joinpath(tempname(), "again_unwrapped.nii")
unwrapping_main([phasefile_me, "-o", phasefile_me_uw, "-t", "[2,4,6]"])
unwrapping_main([phasefile_me_uw, "-o", phasefile_me_uw_wrong, "-t", "[2,4,6]"])
unwrapping_main([phasefile_me_uw, "-o", phasefile_me_uw_again, "-t", "[2,4,6]", "--no-rescale"])

@test readphase(phasefile_me_uw_again; rescale=false).raw == readphase(phasefile_me_uw; rescale=false).raw
@test readphase(phasefile_me_uw_wrong; rescale=false).raw != readphase(phasefile_me_uw; rescale=false).raw

## test ROMEO output files
testpath = joinpath(tmpdir, "test_name_1")
unwrapping_main([phasefile_1eco, "-o", testpath])
@test isfile(joinpath(testpath, "unwrapped.nii"))

testpath = joinpath(tmpdir, "test_name_2")
fn = joinpath(testpath, "unwrap_name.nii")
unwrapping_main([phasefile_1eco, "-o", fn])
@test isfile(fn)

testpath = joinpath(tmpdir, "test_name_2")
gz_fn = joinpath(testpath, "unwrap_name.nii.gz")
unwrapping_main([phasefile_1eco, "-o", gz_fn])
@test isfile(gz_fn)

## test .gz input file
unwrapping_main([gz_fn, "-o", joinpath(testpath, "gz_read_test.nii")])
unwrapping_main([gz_fn, "-m", gz_fn, "-o", joinpath(testpath, "gz_read_test.nii")])


## test mcpc3ds output files
testpath = joinpath(tmpdir, "test5d")
unwrapping_main([phasefile_me_5D, "-o", testpath, "-m", magfile_5D, "-t", "[2,4,6]"])
@test isfile(joinpath(testpath, "combined_mag.nii"))
@test isfile(joinpath(testpath, "combined_phase.nii"))

testpath = joinpath(tmpdir, "test4d")
unwrapping_main([phasefile_me, "-o", testpath, "-m", magfile_me, "-t", "[2,4,6]", "--phase-offset-correction"])
@test isfile(joinpath(testpath, "corrected_phase.nii"))

## test B0 output files
testpath = joinpath(tmpdir, "testB0_1")
unwrapping_main([phasefile_me, "-o", testpath, "-m", magfile_me, "-t", "[2,4,6]", "-B"])
@test isfile(joinpath(testpath, "B0.nii"))

testpath = joinpath(tmpdir, "testB0_2")
name = "B0_output"
unwrapping_main([phasefile_me, "-o", testpath, "-m", magfile_me, "-t", "[2,4,6]", "-B", name])
@test isfile(joinpath(testpath, "$name.nii"))

## TODO add and test homogeneity corrected output

## test quality map
unwrapping_main([phasefile_me, "-m", magfile_me, "-o", tmpdir, "-t", "[2,4,6]", "-qQ"])
fns = joinpath.(tmpdir, ["quality.nii", ("quality_$i.nii" for i in 1:4)...])
for i in 1:length(fns), j in i+1:length(fns)
    @test niread(fns[i]).raw != niread(fns[j]).raw
end

cd(original_path)
GC.gc()
rm(tmpdir, recursive=true)

end

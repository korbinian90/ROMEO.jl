@testset "Threading" begin

# Regression test: `unwrap_individual!` used to share one `Dict` of keyword
# arguments across its `Threads.@threads` loop and write `args[:mag]` into it from
# every iteration, so echo i could be unwrapped using echo j's magnitude. Results
# were non-deterministic on any multi-threaded Julia; on this dataset repeated
# identical calls differed by up to a full 2pi.
#
# Needs real data: on smooth synthetic input the threads finish too quickly to
# interleave and the race does not reproduce. It also needs
# `Threads.nthreads() > 1` - on a single-threaded worker the test still runs but
# cannot fail. Detection is probabilistic, as it is for any race; against the
# unfixed code on 4 threads this caught it in roughly half of the repetitions,
# which is why there are several.

phase = Float64.(niread(joinpath("data", "small", "Phase.nii")).raw)
mag = Float64.(niread(joinpath("data", "small", "Mag.nii")).raw)
TEs = [4.0, 8.0, 12.0]

reference = unwrap_individual(phase; TEs, mag)
for _ in 1:5
    @test unwrap_individual(phase; TEs, mag) == reference
end

end

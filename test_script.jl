include("executable/unwrapping/UnwrappingExecutable.jl")

folder = "C:/ROMEO"
output = joinpath(folder, "output")
phase = joinpath(folder, "ph.nii")
mag = joinpath(folder, "mag.nii")

## the possible arguments can be found in the file executable/unwrapping/UnwrappingExecutable.jl
@time msg = unwrapping_main([phase, "-m", mag, "-o", output, "-t", "[4, 8]"])
println(msg)
GC.gc()

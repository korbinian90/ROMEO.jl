include("UnwrappingExecutable.jl")
# 7T EPI patients
d = "C:/ROMEO"
phase = joinpath(d, "ph.nii")
mag = joinpath(d, "mag.nii")
output = joinpath(d, "output")
@time unwrapping_main([phase, "-o", output, "-m", mag, "-k", "nomask"])
GC.gc()

@testset "voxelquality" begin

phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase4D = niread(phasefile).raw
mag4D = niread(magfile).raw

qm1 = voxelquality(phase4D[:,:,:,1])
qm2 = voxelquality(phase4D[:,:,:,1]; mag=mag4D[:,:,:,1])
qm3 = voxelquality(phase4D; mag4D, TEs=[4,8,12])

@test qm1 != qm2
@test qm1 != qm3
@test qm2 != qm3

end

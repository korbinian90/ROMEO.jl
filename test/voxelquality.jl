@testset "voxelquality" begin

phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase4D = niread(phasefile).raw
mag4D = niread(magfile).raw

qm1 = romeovoxelquality(phase[:,:,:,1])
qm2 = romeovoxelquality(phase[:,:,:,1]; mag=mag[:,:,:,1])
qm3 = romeovoxelquality(phase; mag, TEs=[4,8,12])

@test qm1 != qm2
@test qm1 != qm3
@test qm2 != qm3

end

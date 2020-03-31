using BenchmarkTools

phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
echo = 3
phase = niread(phasefile).raw[:,:,:,echo]
mag = niread(magfile).raw[:,:,:,echo]

big(I, factor=3) = repeat(I; outer=factor .* [1,1,1])
bigphase = big(phase)
bigmag = big(mag)
@btime ROMEO.calculateweights($bigphase; mag=$bigmag)
@btime unwrap($bigphase; mag=$bigmag)

# akku - high power
 # 450ms
 # 1.4s

 # PC
  # 360ms
  # 1.1s

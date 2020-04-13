using BenchmarkTools

phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
echo = 3
phase4D = niread(phasefile).raw
mag4D = niread(magfile).raw
phase = phase4D[:,:,:,echo]
mag = mag4D[:,:,:,echo]

function big(I, factor=3)
    shape = ones(Int, ndims(I))
    shape[1:3] .= factor
    return repeat(I; outer=shape)
end
bigphase = big(phase)
bigmag = big(mag)
@btime ROMEO.calculateweights($bigphase; mag=$bigmag)
@btime unwrap($bigphase; mag=$bigmag)


@btime unwrap($phase; mag=$mag)
@btime unwrap4d!($phase4D; mag=$mag, TEs=1:3)
@btime unwrap4d!($phase4D; mag=$mag4D, TEs=1:3)
@btime unwrap4d!($big(phase4D); mag=$big(mag4D), TEs=1:3)


# akku - high power
 # 450ms
 # 1.4s

# ROMEO Unwrapping
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://korbinian90.github.io/ROMEO.jl/dev)
[![Build Status](https://github.com/korbinian90/ROMEO.jl/workflows/CI/badge.svg)](https://github.com/korbinian90/ROMEO.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/korbinian90/ROMEO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/ROMEO.jl)

Please cite [ROMEO](https://onlinelibrary.wiley.com/doi/10.1002/mrm.28563) if you use it!

## Getting Started

This repository contains ROMEO 3D/4D unwrapping on arrays.  
For MR data in the NIfTI format, a compiled command line tool is available under [ROMEO](https://github.com/korbinian90/ROMEO) (windows and linux executables; does not require a Julia installation) and otherwise, for opening NIfTI files in Julia [NIfTI.jl](https://github.com/JuliaIO/NIfTI.jl) or [MriResearchTools.jl](https://github.com/korbinian90/MriResearchTools.jl) can be helpful.

### Usage

```julia
using ROMEO
unwrapped = unwrap(phasedata3D; mag=magdata3D)
```

or via MriResearchTools:

```julia
using MriResearchTools
phase4D = readphase("Phase.nii") # 4D phase in NIfTI format
unwrapped = unwrap(phase4D; TEs=[4,8,12])
```

## Different Use Cases
### Multi-Echo
If multi-echo data is available, supplying ROMEO with multi-echo information should improve the unwrapping accuracy. The same is true for magnitude information.

### Phase Offsets
If the multi-echo data contains **large phase offsets** (phase at echo time zero), default template unwrapping might fail. Setting the `individual-unwrapping` flag is a solution, as it performs spatial unwrapping for each echo instead. The computed B0 map is not corrected for remaining phase offsets.

For proper handling, the phase offsest can be removed using `mcpc3ds` from `MriResearchTools`. This works for monopolar and bipolar data, already combined or uncombined channels. However, if the phase is already "corrupted" by other coil combination algorithms, it might not be possible to estimate and remove the phase offsets.

### Repeated Measurements (EPI)
4D data with an equal echo time for all volumes should be unwrapped as 4D for best accuracy and temporal stability. The echo times can be set to `TEs=ones(size(phase,4))`

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/ROMEO.jl/blob/master/LICENSE) for details

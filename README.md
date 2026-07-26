# ROMEO Unwrapping

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://korbinian90.github.io/ROMEO.jl/dev)
[![Build Status](https://github.com/korbinian90/ROMEO.jl/workflows/CI/badge.svg)](https://github.com/korbinian90/ROMEO.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/korbinian90/ROMEO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/ROMEO.jl)

Please cite [ROMEO](https://onlinelibrary.wiley.com/doi/10.1002/mrm.28563) if you use it!

## Getting Started

This repository contains ROMEO 3D/4D unwrapping on arrays and a command line interface for MR data in NIfTI format.  
A compiled command line tool is available under [ROMEO](https://github.com/korbinian90/ROMEO) (windows and linux executables; does not require a Julia installation) and otherwise, for opening NIfTI files in Julia [NIfTI.jl](https://github.com/JuliaIO/NIfTI.jl) or [MriResearchTools.jl](https://github.com/korbinian90/MriResearchTools.jl) can be helpful.

### Usage - command line

Install Julia 1.9 or newer (https://julialang.org)  
Copy the file [romeo.jl](https://github.com/korbinian90/ROMEO.jl/blob/master/romeo.jl) from this repository to a convenient location. An alias for `romeo` as `julia <path-to-file>/romeo.jl` might be useful.

```bash
    $ julia <path-to-file>/romeo.jl phase.nii -m mag.nii -t [2.1,4.2,6.3] -o results
```

On the first run, the dependencies will be installed automatically.

For an extended explanation of the command line interface see [ROMEO](https://github.com/korbinian90/ROMEO)  

### Usage - Julia

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

### Function Reference
https://korbinian90.github.io/ROMEO.jl/dev

## Different Use Cases
### Multi-Echo
If multi-echo data is available, supplying ROMEO with multi-echo information should improve the unwrapping accuracy. The same is true for magnitude information.

### Phase Offsets
If the multi-echo data contains **large phase offsets** (phase at echo time zero), default template unwrapping might fail. Setting the `individual-unwrapping` flag is a solution, as it performs spatial unwrapping for each echo instead. The computed B0 map is not corrected for remaining phase offsets.

For proper handling, the phase offsest can be removed using `mcpc3ds` from `MriResearchTools`. This works for monopolar and bipolar data, already combined or uncombined channels. However, if the phase is already "corrupted" by other coil combination algorithms, it might not be possible to estimate and remove the phase offsets.

### Repeated Measurements (EPI)
4D data with an equal echo time for all volumes should be unwrapped as 4D for best accuracy and temporal stability. The echo times can be set to `TEs=ones(size(phase,4))`

### Phase Singularities around Veins (experimental)
In high resolution data (e.g. 7T, 0.3-0.6 mm, TE >= 20 ms) the complex signal can cross zero at a vessel wall, where two compartments with roughly opposite phase cancel each other. At such a point the wrapped phase has a non-zero circulation (a residue), which is a topological obstruction: no congruent unwrapping exists. In 3D the residue field is divergence free, so the residues form closed lines, around a vein typically rings that encircle the vessel. ROMEO stays congruent and therefore has to keep the resulting branch cut, it only moves it to the worst weights.

Detection is parameter free and can be run on the wrapped input or on the unwrapped output:

```julia
using ROMEO
unwrapped = unwrap(phase; mag)
s = detect_singularities(unwrapped)
mask = singularity_mask(s)   # voxels on a vortex line
cuts = branchcuts(unwrapped) # faces with a remaining 2π jump
```

The singularity mask alone is already useful as a data term weight for a QSM inversion. The phase can also be corrected locally, which keeps the result bit-identical outside of small patches around the vortex lines:

```julia
info = fix_singularities!(unwrapped; mag, mode=:cascade)
info.nfixed, info.nskipped, info.nchanged
info.delta # the applied correction, φ_corrected - φ_ROMEO
```

On the command line, `--fix-singularities [off|mask|cascade|smooth|lsq|inpaint]` (default `off`) writes `singularities.nii` and `branchcuts.nii`, and for the correcting modes also `phase_correction.nii`. Loops that are too long, too large or that are not in a signal void are only marked, never modified, which protects real pathology (microbleeds, implants).

### Setting the Template Echo
In certain cases, the phase of the first echo/time-point looks differently than the rest of the acquisition, which can occur due to flow compensation of only the first echo or not having reached the steady state in fMRI. This might cause template unwrapping to fail, as the first echo is chosen as the template by default.  
With the optional argument `template=2`, this can be changed to the second (or any other) echo/time-point.

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/ROMEO.jl/blob/master/LICENSE) for details

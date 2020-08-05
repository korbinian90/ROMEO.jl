# ROMEO Unwrapping
[![Build Status](https://travis-ci.com/korbinian90/ROMEO.jl.svg?branch=master)](https://travis-ci.com/korbinian90/ROMEO.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/korbinian90/ROMEO.jl?svg=true)](https://ci.appveyor.com/project/korbinian90/ROMEO-jl)
[![Codecov](https://codecov.io/gh/korbinian90/ROMEO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/ROMEO.jl)
[![Coveralls](https://coveralls.io/repos/github/korbinian90/ROMEO.jl/badge.svg?branch=master)](https://coveralls.io/github/korbinian90/ROMEO.jl?branch=master)

Please cite [ROMEO bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.24.214551v1.abstract) if you use it! The link will update to the peer reviewed version after it is published.

## Getting Started

This repository contains ROMEO 3D/4D unwrapping on arrays.

For MR data in the NIfTI format, a compiled command line tool is available under [ROMEO](https://github.com/korbinian90/ROMEO) (windows and linux binaries; does not require a Julia installation) and otherwise, for opening NIfTI files in Julia [NIfTI.jl](https://github.com/JuliaIO/NIfTI.jl) or [MriResearchTools.jl](https://github.com/korbinian90/MriResearchTools.jl) can be helpful.

### Usage

```julia
using ROMEO
unwrapped = unwrap(phasedata3D; mag=magdata3D)
```

or via MriResearchTools:

```julia
using MriResearchTools
phase4D = readphase("Phase.nii") # 4D phase in NIfTI format
unwrapped = unwrap(phase4D)
```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/ROMEO.jl/blob/master/LICENSE) for details

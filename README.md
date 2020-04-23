# ROMEO Unwrapping
[![Build Status](https://travis-ci.com/korbinian90/ROMEO.jl.svg?branch=master)](https://travis-ci.com/korbinian90/ROMEO.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/korbinian90/ROMEO.jl?svg=true)](https://ci.appveyor.com/project/korbinian90/ROMEO-jl)
[![Codecov](https://codecov.io/gh/korbinian90/ROMEO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/ROMEO.jl)
[![Coveralls](https://coveralls.io/repos/github/korbinian90/ROMEO.jl/badge.svg?branch=master)](https://coveralls.io/github/korbinian90/ROMEO.jl?branch=master)

Please cite ROMEO if you are applying it in your method.

## Getting Started

This repository contains the ROMEO algorithm for 3D unwrapping.

For the specialized version for MR data in the nifti format (also 4D unwrapping), see [MriResearchTools](https://github.com/korbinian90/MriResearchTools.jl).

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

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/ROMEO.jl/blob/master/LICENSE) for details

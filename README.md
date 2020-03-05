# ROMEO Unwrapping
[![Build Status](https://travis-ci.com/korbinian90/ROMEO.jl.svg?branch=master)](https://travis-ci.com/korbinian90/ROMEO.jl)
[![Codecov](https://codecov.io/gh/korbinian90/ROMEO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/ROMEO.jl)
[![Coveralls](https://coveralls.io/repos/github/korbinian90/ROMEO.jl/badge.svg?branch=master)](https://coveralls.io/github/korbinian90/ROMEO.jl?branch=master)

The version of the code used for the publication of ROMEO [doi-inserted-after-publication]. Please cite ROMEO if you are applying it in your method.

## Getting Started

This repository contains the ROMEO algorithm for 3D unwrapping on arrays.

For the specialized version for MR data in the nifti format (also 4D unwrapping), see https://github.com/korbinian90/MriResearchTools.jl.

### Installing the Julia code
A Julia installation v1.1 or higher is required.

Open the julia REPL inside the the ROMEO.jl folder and type
```julia
julia> ] # enter julia package manager
(v1.x) pkg> add https://github.com/korbinian90/ROMEO.jl
(ROMEO.jl) pkg> # type backspace to get back to the julia REPL
julia>
```

### Usage

```julia
using ROMEO
unwrapped = unwrap(phasedata3D; mag=magdata3D)
```

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/ROMEO.jl/blob/master/LICENSE) for details

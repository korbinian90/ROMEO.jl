# ROMEO Unwrapping

The version of the code used for the publication of ROMEO [doi-inserted-after-publication]. Please cite ROMEO if you are applying it in your method.

## Getting Started

### Prerequisites

Magnitude and Phase images in NIfTI fileformat (4D images with echoes in the 4th dimension).

The datasets used in the publication can be obtained from https://dataverse.harvard.edu/dataverse/ROMEO.

### Compiled Version
Compiled versions for windows and linux are attached to the release v1.0 in github.

Example usage:
```shell
romeo ph.nii -m mag.ii -k nomask -o outputdir
```

Example usage of the compiled version from MATLAB:
```Matlab
!romeo ph.nii -m mag.ii -k nomask -o outputdir
```

### Installing the Julia code
A Julia installation v1.1 or higher is required.

Clone the git repository and checkout the publication branch.
Alternatively, download the publication branch as a zip file end extract.

```shell
cd <folder-to-download-romeo>
git clone https://github.com/korbinian90/ROMEO.jl
cd ROMEO.jl
git checkout publication
```

Open the julia REPL inside the the ROMEO.jl folder and type
```julia
julia> ] # enter julia package manager
(v1.x) pkg> activate . # Activate the ROMEO.jl environment
(ROMEO.jl) pkg> instantiate # Install all dependencies
(ROMEO.jl) pkg> # type backspace to get back to the julia REPL
julia>
```

### Usage
ROMEO can be used via the julia REPL, via julia from the command line or via a compiled program from the command line.

Julia REPL usage:

Set the parameters in the file test_script.jl run in the julia REPL
```julia
include("test_script.jl")
```

Julia command line usage:
```shell
julia test.jl ph.nii -m mag.nii -o outputdir
```

Compiling command line program:

Set the buildpath in executable/compile.jl and run in the julia REPL
```julia
include("executable/compile.jl")
```

The executable is called “romeo” on linux and “romeo.exe” on windows.

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/ROMEO.jl/blob/publication/LICENSE) for details

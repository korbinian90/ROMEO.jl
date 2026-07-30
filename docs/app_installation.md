# Installing ROMEO as a Julia App

ROMEO can be installed as a Julia Pkg app, which provides convenient command-line access to the `unwrapping_main` function.

## Installation

### Using Julia 1.9 or later:

```julia
using Pkg
Pkg.add("ROMEO", app=true)
```

Or from the Pkg REPL (press `]` in Julia):

```
pkg> app add ROMEO
```

This will install the `romeo` command-line tool to `~/.julia/bin/romeo`. Make sure `~/.julia/bin` is in your system PATH for easy access.

## Usage

After installation, you can run ROMEO from the command line:

```bash
romeo --help
romeo --phase phase.nii.gz --magnitude mag.nii.gz --output unwrapped.nii.gz
```

## Updating

To update the app:

```julia
using Pkg
Pkg.update("ROMEO")
```

Or:

```
pkg> app update ROMEO
```

## Uninstallation

To remove the app:

```julia
using Pkg
Pkg.rm("ROMEO", app=true)
```

Or:

```
pkg> app rm ROMEO
```

## Alternative: Direct Script Usage

If you prefer not to use the app installation, you can still use the provided `romeo.jl` script as described in the main README.

function getargs(args::AbstractVector, version)
    if isempty(args)
        args = ["--help"]
    else
        if !('-' in args[1]) prepend!(args, Ref("-p")) end # if phase is first without -p
        if length(args) >= 2 && !("-p" in args || "--phase" in args) && !('-' in args[end-1]) # if phase is last without -p
            insert!(args, length(args), "-p")
        end
    end
    s = ArgParseSettings(;
        exc_handler=exception_handler,
        add_version=true,
        version,
        )
    @add_arg_table! s begin
        "--phase", "-p"
            help = "The phase image that should be unwrapped"
        "--magnitude", "-m"
            help = "The magnitude image (better unwrapping if specified)"
        "--output", "-o"
            help = "The output path or filename"
            default = "unwrapped.nii"
        "--echo-times", "-t"
            help = """The echo times in [ms] required for temporal unwrapping 
                specified in array or range syntax (eg. "[1.5,3.0]" or 
                "3.5:3.5:14"). For identical echo times, "-t epi" can be
                used with the possibility to specify the echo time as
                e.g. "-t epi 5.3" (for B0 calculation)."""
            nargs = '+'
        "--mask", "-k"
            help = """nomask | qualitymask <threshold> | robustmask | <mask_file>.
                <threshold>=0.1 for qualitymask in [0;1]"""
            default = ["robustmask"]
            nargs = '+'
        "--mask-unwrapped", "-u"
            help = """Apply the mask on the unwrapped result. If mask is 
                "nomask", sets it to "robustmask"."""
            action = :store_true
        "--unwrap-echoes", "-e"
            help = "Load only the specified echoes from disk"
            default = [":"]
            nargs = '+'
        "--weights", "-w"
            help = """romeo | romeo2 | romeo3 | romeo4 | romeo6 |
                bestpath | <4d-weights-file> | <flags>.
                <flags> are up to 6 bits to activate individual weights
                (eg. "1010"). The weights are (1)phasecoherence
                (2)phasegradientcoherence (3)phaselinearity (4)magcoherence
                (5)magweight (6)magweight2"""
            default = "romeo"
        "--compute-B0", "-B"
            help = """Calculate combined B0 map in [Hz].
                Supports the B0 output filename as optional input.
                This activates MCPC3Ds phase offset correction (monopolar)
                for multi-echo data."""
            default = ""
            nargs = '?'
            constant = "B0"
        "--phase-offset-correction"
            help = """on | off | bipolar.
                Applies the MCPC3Ds method to perform phase offset
                determination and removal (for multi-echo).
                "bipolar" removes eddy current artefacts (requires >= 3 echoes)."""
            default = "default off"
            nargs = '?'
            constant = "on"
        "--phase-offset-smoothing-sigma-mm"
            help = """default: [7,7,7]
                Only applied if phase-offset-correction is activated. The given
                sigma size is divided by the voxel size from the nifti phase
                file to obtain a smoothing size in voxels. A value of [0,0,0]
                deactivates phase offset smoothing (not recommended)."""
            nargs = '+'
        "--write-phase-offsets"
            help = "Saves the estimated phase offsets to the output folder"
            action = :store_true
        "--individual-unwrapping", "-i"
            help = """Unwraps the echoes individually (not temporal).
                This might be necessary if there is large movement
                (timeseries) or phase-offset-correction is not
                applicable."""
            action = :store_true
        "--template"
            help = """Template echo that is spatially unwrapped and used for
                temporal unwrapping"""
            arg_type = Int
            default = 1
        "--no-mmap", "-N"
            help = """Deactivate memory mapping. Memory mapping might cause
                problems on network storage"""
            action = :store_true
        "--no-rescale"
            help = """Deactivate rescaling of input images. By default the
                input phase is rescaled to the range [-π;π]. This option
                allows inputting already unwrapped phase images without
                manually wrapping them first."""
            action = :store_true
        "--threshold"
            help = """<maximum number of wraps>.
                Threshold the unwrapped phase to the maximum number of wraps
                and sets exceeding values to 0"""
            arg_type = Float64
            default = Inf
        "--verbose", "-v"
            help = "verbose output messages"
            action = :store_true
        "--correct-global", "-g"
            help = """Phase is corrected to remove global n2π phase offset. The
                median of phase values (inside mask if given) is used to
                calculate the correction term"""
            action = :store_true
        "--write-quality", "-q"
            help = """Writes out the ROMEO quality map as a 3D image with one
                value per voxel"""
            action = :store_true
        "--write-quality-all", "-Q"
            help = """Writes out an individual quality map for each of the
                ROMEO weights."""
            action = :store_true
        "--max-seeds", "-s"
            help = """EXPERIMENTAL! Sets the maximum number of seeds for
                unwrapping. Higher values allow more seperated regions."""
            arg_type = Int
            default = 1
        "--merge-regions"
            help = """EXPERIMENTAL! Spatially merges neighboring regions after
                unwrapping."""
            action = :store_true
        "--correct-regions"
            help = """EXPERIMENTAL! Performed after merging. Brings the median
                of each region closest to 0 (mod 2π)."""
            action = :store_true
        "--wrap-addition"
            help = """[0;π] EXPERIMENTAL! Usually the true phase difference of
                neighboring voxels cannot exceed π to be able to unwrap
                them. This setting increases the limit and uses 'linear
                unwrapping' of 3 voxels in a line. Neighbors can have
                (π + wrap-addition) phase difference."""
            arg_type = Float64
            default = 0.0
        "--temporal-uncertain-unwrapping"
            help = """EXPERIMENTAL! Uses spatial unwrapping on voxels that have
                high uncertainty values after temporal unwrapping."""
            action = :store_true

    end
    return parse_args(args, s)
end

function exception_handler(settings::ArgParseSettings, err, err_code::Int=1)
    if err == ArgParseError("too many arguments")
        println(stderr,
            """wrong argument formatting!"""
        )
    end
    ArgParse.default_handler(settings, err, err_code)
end

function getechoes(settings, neco)
    echoes = eval(Meta.parse(join(settings["unwrap-echoes"], " ")))
    if echoes isa Int
        echoes = [echoes]
    elseif echoes isa Matrix
        echoes = echoes[:]
    end
    echoes = (1:neco)[echoes] # expands ":"
    if (length(echoes) == 1) echoes = echoes[1] end
    return echoes
end

function getTEs(settings, neco, echoes)
    if isempty(settings["echo-times"])
        if neco == 1 || length(echoes) == 1
            return 1
        else
            error("multi-echo data is used, but no echo times are given. Please specify the echo times using the -t option.")
        end
    end
    TEs = if settings["echo-times"][1] == "epi"
        ones(neco) .* if length(settings["echo-times"]) > 1; parse(Float64, settings["echo-times"][2]) else 1 end
    else
        eval(Meta.parse(join(settings["echo-times"], " ")))
    end
    if TEs isa AbstractMatrix
        TEs = TEs[:]
    end
    if 1 < length(TEs) == neco
        TEs = TEs[echoes]
    end
    return TEs
end

function get_phase_offset_smoothing_sigma(settings)
    if isempty(settings["phase-offset-smoothing-sigma-mm"])
        return (7,7,7)
    end
    return eval(Meta.parse(join(settings["phase-offset-smoothing-sigma-mm"], " ")))[:]
end

function parseweights(settings)
    if isfile(settings["weights"]) && splitext(settings["weights"])[2] != ""
        return UInt8.(niread(settings["weights"]))
    else
        try
            reform = "Bool[$(join(collect(settings["weights"]), ','))]"
            flags = falses(6)
            flags_tmp = eval(Meta.parse(reform))
            flags[1:length(flags_tmp)] = flags_tmp 
            return flags
        catch
            return Symbol(settings["weights"])
        end
    end
end

function saveconfiguration(writedir, settings, args, version)
    writedir = abspath(writedir)
    open(joinpath(writedir, "settings_romeo.txt"), "w") do io
        for (fname, val) in settings
            if !(val isa AbstractArray || fname == "header")
                println(io, "$fname: " * string(val))
            end
        end
        println(io, """Arguments: $(join(args, " "))""")
        println(io, "RomeoApp version: $version")
    end
end

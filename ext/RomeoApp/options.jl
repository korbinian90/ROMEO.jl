const CLI = ROMEO.CLI

const OPTIONS = [
    CLI.Option("--phase", "-p", "The phase image that should be unwrapped"),
    CLI.Option("--magnitude", "-m", "The magnitude image (better unwrapping if specified)"),
    CLI.Option("--output", "-o", """The output path or filename. The unwrapped phase is
        written in radians, keeping the wraps that were removed as
        multiples of 2pi, so its range is wider than the input.
        (default: unwrapped.nii)"""),
    CLI.Option("--echo-times", "-t", """The echo times in [ms] required for temporal unwrapping
        specified in array or range syntax (eg. "[1.5,3.0]" or
        "3.5:3.5:14"). For identical echo times, "-t epi" can be
        used with the possibility to specify the echo time as
        e.g. "-t epi 5.3" (for B0 calculation)."""; nargs=:many),
    CLI.Option("--mask", "-k", """nomask | qualitymask <threshold> | robustmask | <mask_file>.
        <threshold>=0.1 for qualitymask in [0;1] (default: robustmask)"""; nargs=:many),
    CLI.Option("--mask-unwrapped", "-u", """Apply the mask on the unwrapped result. If mask is
        "nomask", sets it to "robustmask"."""; nargs=:none),
    CLI.Option("--unwrap-echoes", "-e", "Load only the specified echoes from disk (default: :)"; nargs=:many),
    CLI.Option("--weights", "-w", """romeo | romeo2 | romeo3 | romeo4 | romeo6 |
        bestpath | <4d-weights-file> | <flags>.
        <flags> are up to 6 bits to activate individual weights
        (eg. "1010"). The weights are (1)phasecoherence
        (2)phasegradientcoherence (3)phaselinearity (4)magcoherence
        (5)magweight (6)magweight2 (default: romeo)"""),
    CLI.Option("--compute-B0", "-B", """Calculate combined B0 map in [Hz].
        Supports the B0 output filename as optional input.
        This activates MCPC3Ds phase offset correction (monopolar)
        for multi-echo data."""; nargs=:optional),
    CLI.Option("--B0-phase-weighting", "", """phase_snr | phase_var | average | TEs | mag | simulated_mag
        Set the weighting for the B0 calculation. (default: phase_snr)"""),
    CLI.Option("--phase-offset-correction", "", """on | off | bipolar.
        Applies the MCPC3Ds method to perform phase offset
        determination and removal (for multi-echo).
        "bipolar" removes eddy current artefacts (requires >= 3 echoes)."""; nargs=:optional),
    CLI.Option("--phase-offset-smoothing-sigma-mm", "", """default: [7,7,7]
        Only applied if phase-offset-correction is activated. The given
        sigma size is divided by the voxel size from the nifti phase
        file to obtain a smoothing size in voxels. A value of [0,0,0]
        deactivates phase offset smoothing (not recommended).
        The 7mm default assumes a whole-head field of view. On a small
        field of view (ex vivo, small animal, a single slab) 7mm covers
        a large part of the object, which smooths away the phase offset
        structure it is meant to capture and leaves a gradient across
        the fieldmap. Reduce it to a few voxels in that case, and use 0
        for a dimension with a single slice."""; nargs=:many),
    CLI.Option("--write-phase-offsets", "", "Saves the estimated phase offsets to the output folder"; nargs=:none),
    CLI.Option("--individual-unwrapping", "-i", """Unwraps the echoes individually (not temporal).
        This might be necessary if there is large movement
        (timeseries) or phase-offset-correction is not
        applicable."""; nargs=:none),
    CLI.Option("--template", "", """Template echo that is spatially unwrapped and used for
        temporal unwrapping (default: 1)"""),
    CLI.Option("--no-mmap", "-N", """Deactivate memory mapping. Memory mapping might cause
        problems on network storage"""; nargs=:none),
    CLI.Option("--no-phase-rescale", "", """Deactivate rescaling of input images. By default the
        input phase is rescaled to the range [-π;π]. This option
        allows inputting already unwrapped phase images without
        manually wrapping them first."""; nargs=:none),
    CLI.Option("--no-rescale", "", "Same as --no-phase-rescale"; nargs=:none),
    CLI.Option("--fix-ge-phase", "", """GE systems write corrupted phase output (slice jumps).
        This option fixes the phase problems."""; nargs=:none),
    CLI.Option("--threshold", "", """<maximum number of wraps>.
        Threshold the unwrapped phase to the maximum number of wraps
        and sets exceeding values to 0 (default: Inf)"""),
    CLI.Option("--verbose", "-v", "verbose output messages"; nargs=:none),
    CLI.Option("--correct-global", "-g", """Phase is corrected to remove global n2π phase offset. The
        median of phase values (inside mask if given) is used to
        calculate the correction term. This also corrects multi-echo
        phase for individual unwrapping, and might require MCPC3Ds
        phase offset correction."""; nargs=:none),
    CLI.Option("--write-quality", "-q", """Writes out the ROMEO quality map as a 3D image with one
        value per voxel"""; nargs=:none),
    CLI.Option("--write-quality-all", "-Q", """Writes out an individual quality map for each of the
        ROMEO weights."""; nargs=:none),
    CLI.Option("--max-seeds", "-s", """EXPERIMENTAL! Sets the maximum number of seeds for
        unwrapping. Higher values allow more seperated regions. (default: 1)"""),
    CLI.Option("--merge-regions", "", """EXPERIMENTAL! Spatially merges neighboring regions after
        unwrapping."""; nargs=:none),
    CLI.Option("--correct-regions", "", """EXPERIMENTAL! Performed after merging. Brings the median
        of each region closest to 0 (mod 2π)."""; nargs=:none),
    CLI.Option("--wrap-addition", "", """[0;π] EXPERIMENTAL! Usually the true phase difference of
        neighboring voxels cannot exceed π to be able to unwrap
        them. This setting increases the limit and uses 'linear
        unwrapping' of 3 voxels in a line. Neighbors can have
        (π + wrap-addition) phase difference. (default: 0.0)"""),
    CLI.Option("--temporal-uncertain-unwrapping", "", """EXPERIMENTAL! Uses spatial unwrapping on voxels that have
        low quality values after temporal unwrapping. A higher threshold
        leads to more voxels being spatially unwrapped. The range of the
        threshold is [0;1] (default: 0, 0.5 when given without value)"""; nargs=:optional),
]

# The command line, resolved. Fields after `temporal_uncertain_unwrapping` are
# filled in while running and end up in the settings record.
Base.@kwdef mutable struct RomeoOptions
    phase::String = ""
    magnitude::String = ""
    output::String = "unwrapped.nii"
    echo_times::Vector{String} = String[]
    mask::Vector{String} = ["robustmask"]
    mask_unwrapped::Bool = false
    unwrap_echoes::Vector{String} = [":"]
    weights::String = "romeo"
    compute_B0::String = ""
    B0_phase_weighting::String = "phase_snr"
    phase_offset_correction::String = "default off"
    phase_offset_smoothing_sigma_mm::Vector{String} = String[]
    write_phase_offsets::Bool = false
    individual_unwrapping::Bool = false
    template::Int = 1
    no_mmap::Bool = false
    no_phase_rescale::Bool = false
    fix_ge_phase::Bool = false
    threshold::Float64 = Inf
    verbose::Bool = false
    correct_global::Bool = false
    write_quality::Bool = false
    write_quality_all::Bool = false
    max_seeds::Int = 1
    merge_regions::Bool = false
    correct_regions::Bool = false
    wrap_addition::Float64 = 0.0
    temporal_uncertain_unwrapping::Float64 = 0.0
    filename::String = "unwrapped"
    number_of_echoes::Int = 0
    multi_channel::Bool = false
    echoes::Vector{Int} = Int[]
    TEs::Union{Vector{Int},Vector{Float64}} = Int[]
end

one(v::CLI.Values, key, default::String) = haskey(v, key) ? v[key][1] : default
many(v::CLI.Values, key, default::Vector{String}) = haskey(v, key) ? v[key] : default
flag(v::CLI.Values, key) = haskey(v, key)
optional(v::CLI.Values, key, default::String, constant::String) = haskey(v, key) ? (isempty(v[key]) ? constant : v[key][1]) : default

function RomeoOptions(v::CLI.Values)
    tuu = optional(v, "temporal-uncertain-unwrapping", "0", "0.5")
    return RomeoOptions(;
        phase = one(v, "phase", ""),
        magnitude = one(v, "magnitude", ""),
        output = one(v, "output", "unwrapped.nii"),
        echo_times = many(v, "echo-times", String[]),
        mask = many(v, "mask", ["robustmask"]),
        mask_unwrapped = flag(v, "mask-unwrapped"),
        unwrap_echoes = many(v, "unwrap-echoes", [":"]),
        weights = one(v, "weights", "romeo"),
        compute_B0 = optional(v, "compute-B0", "", "B0"),
        B0_phase_weighting = one(v, "B0-phase-weighting", "phase_snr"),
        phase_offset_correction = optional(v, "phase-offset-correction", "default off", "on"),
        phase_offset_smoothing_sigma_mm = many(v, "phase-offset-smoothing-sigma-mm", String[]),
        write_phase_offsets = flag(v, "write-phase-offsets"),
        individual_unwrapping = flag(v, "individual-unwrapping"),
        template = parse(Int, one(v, "template", "1")),
        no_mmap = flag(v, "no-mmap"),
        no_phase_rescale = flag(v, "no-phase-rescale") || flag(v, "no-rescale"),
        fix_ge_phase = flag(v, "fix-ge-phase"),
        threshold = parse(Float64, one(v, "threshold", "Inf")),
        verbose = flag(v, "verbose"),
        correct_global = flag(v, "correct-global"),
        write_quality = flag(v, "write-quality"),
        write_quality_all = flag(v, "write-quality-all"),
        max_seeds = parse(Int, one(v, "max-seeds", "1")),
        merge_regions = flag(v, "merge-regions"),
        correct_regions = flag(v, "correct-regions"),
        wrap_addition = parse(Float64, one(v, "wrap-addition", "0.0")),
        temporal_uncertain_unwrapping = parse(Float64, tuu),
    )
end

# The settings record, with the names of the command line options
function settings_record(o::RomeoOptions)
    f = CLI.format
    return Dict{String,String}(
        "phase" => f(o.phase),
        "magnitude" => f(o.magnitude),
        "output" => f(o.output),
        "echo-times" => f(o.echo_times),
        "mask" => f(o.mask),
        "mask-unwrapped" => f(o.mask_unwrapped),
        "unwrap-echoes" => f(o.unwrap_echoes),
        "weights" => f(o.weights),
        "compute-B0" => f(o.compute_B0),
        "B0-phase-weighting" => f(o.B0_phase_weighting),
        "phase-offset-correction" => f(o.phase_offset_correction),
        "phase-offset-smoothing-sigma-mm" => f(o.phase_offset_smoothing_sigma_mm),
        "write-phase-offsets" => f(o.write_phase_offsets),
        "individual-unwrapping" => f(o.individual_unwrapping),
        "template" => f(o.template),
        "no-mmap" => f(o.no_mmap),
        "no-phase-rescale" => f(o.no_phase_rescale),
        "fix-ge-phase" => f(o.fix_ge_phase),
        "threshold" => f(o.threshold),
        "verbose" => f(o.verbose),
        "correct-global" => f(o.correct_global),
        "write-quality" => f(o.write_quality),
        "write-quality-all" => f(o.write_quality_all),
        "max-seeds" => f(o.max_seeds),
        "merge-regions" => f(o.merge_regions),
        "correct-regions" => f(o.correct_regions),
        "wrap-addition" => f(o.wrap_addition),
        "temporal-uncertain-unwrapping" => f(o.temporal_uncertain_unwrapping),
        "filename" => f(o.filename),
        "number-of-echoes" => f(o.number_of_echoes),
        "multi-channel" => f(o.multi_channel),
        "echoes" => format_scalar_or_vector(o.echoes),
        "TEs" => format_scalar_or_vector(o.TEs),
    )
end

# a single echo or echo time is recorded as the number itself
format_scalar_or_vector(v::AbstractVector) = length(v) == 1 ? CLI.format(v[1]) : CLI.format(v)

function getargs(args::AbstractVector, version)
    if isempty(args)
        args = ["--help"]
    else
        # `startswith`, not `'-' in`: the latter asked whether a hyphen appears
        # anywhere in the string, so any filename containing one looked like a
        # flag and did not get `-p` prepended. BIDS names are made of hyphens
        # (sub-01_part-phase_MEGRE.nii), so the documented form
        # `romeo phase.nii -t ...` failed for most real datasets.
        if !startswith(args[1], '-') prepend!(args, Ref("-p")) end # if phase is first without -p
        if length(args) >= 2 && !("-p" in args || "--phase" in args) && !startswith(args[end-1], '-') # if phase is last without -p
            insert!(args, length(args), "-p")
        end
    end
    spec = CLI.Spec("romeo", string(version), OPTIONS)
    values = CLI.parse(spec, args)
    values === nothing && return nothing
    return RomeoOptions(values)
end

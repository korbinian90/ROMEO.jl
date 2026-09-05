function ROMEO.unwrapping_main(args; version=package_version(ROMEO))
    opts = try
        getargs(args, version)
    catch e
        e isa ArgumentError || rethrow()
        msg = e.msg
        msg isa String || (msg = "wrong arguments")
        CLI.print_stderr("romeo: " * msg * "\nwrong argument formatting! See romeo --help\n")
        return 1
    end
    opts === nothing && return 0 # --help or --version
    return run_romeo(opts, args, version)
end

function run_romeo(opts::RomeoOptions, args, version)
    phase, mag, hdr = load_data_and_resolve_args!(opts)

    mkpath(opts.output)
    saveconfiguration(opts, args, version)

    ## Perform phase offset correction with MCPC-3D-S (and in case of 5D coil combination)
    phase, mag = phase_offset_correction!(opts, phase, mag, hdr)

    # A single echo is unwrapped as 3D data, several echoes as 4D data.
    if length(opts.echoes) == 1
        unwrap_and_write(opts, select_echo(phase, opts.echoes[1]), select_echo(mag, opts.echoes[1]), hdr, opts.TEs)
    else
        unwrap_and_write(opts, select_echoes(phase, opts.echoes), select_echoes(mag, opts.echoes), hdr, opts.TEs)
    end
    return 0
end

size3(a) = (size(a, 1), size(a, 2), size(a, 3))
select_echo(a::AbstractArray{<:Any,4}, i::Int) = a[:,:,:,i]
select_echo(::Nothing, i) = nothing
select_echoes(a::AbstractArray{<:Any,4}, echoes::Vector{Int}) = a[:,:,:,echoes]
select_echoes(::Nothing, echoes) = nothing
as4d(a::AbstractArray{<:Any,5}) = reshape(a, Base.front(size(a)))
as4d(::Nothing) = nothing

function load_data_and_resolve_args!(opts::RomeoOptions)
    opts.filename = "unwrapped"
    if endswith(opts.output, ".nii") || endswith(opts.output, ".nii.gz")
        opts.filename = basename(opts.output)
        opts.output = dirname(opts.output)
    end

    if opts.weights == "romeo"
        opts.weights = isempty(opts.magnitude) ? "romeo4" : "romeo3"
    end

    if opts.mask_unwrapped && opts.mask[1] == "nomask"
        opts.mask[1] = "robustmask"
    end

    if opts.mask[1] == "robustmask" && isempty(opts.magnitude)
        opts.mask[1] = "nomask"
        CLI.warn("robustmask was chosen but no magnitude is available. No mask is used!")
    end

    isempty(opts.phase) && error("no phase image given (-p)")
    phase, hdr = loadphase(opts.phase; rescale=!opts.no_phase_rescale, fix_ge=opts.fix_ge_phase)
    opts.verbose && CLI.info("Phase loaded!")
    mag = isempty(opts.magnitude) ? nothing : first(loadmag(opts.magnitude))
    opts.verbose && mag !== nothing && CLI.info("Mag loaded!")

    neco = opts.number_of_echoes = size(phase, 4)

    # activate phase-offset-correction as default (monopolar)
    opts.multi_channel = size(phase, 5) > 1
    if (!isempty(opts.compute_B0) || opts.multi_channel || opts.phase_offset_correction == "on") && !(opts.phase_offset_correction in ("bipolar", "off"))
        opts.phase_offset_correction = "monopolar"
        opts.verbose && CLI.info("Phase offset correction with MCPC-3D-S set to monopolar")
    end
    if neco == 1
        opts.phase_offset_correction = "off"
        opts.verbose && CLI.info("Phase offset correction with MCPC-3D-S turned off (only one echo)")
    end
    if opts.phase_offset_correction == "default off"
        opts.phase_offset_correction = "off"
    end

    ## Echoes for unwrapping
    opts.echoes = try
        getechoes(opts, neco)
    catch y
        if isa(y, BoundsError)
            error("echoes=$(join(opts.unwrap_echoes, " ")): specified echo out of range! Number of echoes is $neco")
        else
            error("echoes=$(join(opts.unwrap_echoes, " ")) wrongly formatted!")
        end
    end
    opts.verbose && CLI.info("Echoes are " * format_scalar_or_vector(opts.echoes))

    opts.TEs = getTEs(opts, neco, opts.echoes)
    opts.verbose && CLI.info("TEs are " * format_scalar_or_vector(opts.TEs))

    if 1 < length(opts.echoes) && length(opts.echoes) != length(opts.TEs)
        error("Number of chosen echoes is $(length(opts.echoes)) ($neco in .nii data), but $(length(opts.TEs)) TEs were specified!")
    end

    if mag !== nothing && (size3(mag) != size3(phase) || size(mag, 4) < maximum(opts.echoes))
        error("size of magnitude and phase does not match!")
    end

    equal_echo_time = length(opts.TEs) >= 2 && opts.TEs[1] == opts.TEs[2]
    if opts.phase_offset_correction != "off" && equal_echo_time
        CLI.warn("The echo times 1 and 2 ($(CLI.format(opts.TEs))) need to be different for MCPC-3D-S phase offset correction! No phase offset correction performed")
        opts.phase_offset_correction = "off"
    end

    return phase, mag, hdr
end

function getechoes(opts::RomeoOptions, neco)
    parsed = ROMEO.parse_array(opts.unwrap_echoes)
    echoes = if parsed isa Colon
        collect(1:neco)
    elseif parsed isa Int
        [parsed]
    elseif parsed isa Vector{Int}
        parsed
    else
        throw(ArgumentError("echoes must be integers"))
    end
    return (1:neco)[echoes] # throws BoundsError for an echo out of range
end

# A number or a vector of numbers, as a vector
function parse_numbers(strs)
    parsed = ROMEO.parse_array(strs)
    parsed isa Int && return [parsed]
    parsed isa Float64 && return [parsed]
    parsed isa Vector{Int} && return parsed
    parsed isa Vector{Float64} && return parsed
    throw(ArgumentError("\"$(join(strs, " "))\" is not a number or an array of numbers"))
end

function getTEs(opts::RomeoOptions, neco, echoes)
    if isempty(opts.echo_times)
        if neco == 1 || length(echoes) == 1
            return [1]
        else
            error("multi-echo data is used, but no echo times are given. Please specify the echo times using the -t option.")
        end
    end
    TEs = if opts.echo_times[1] == "epi"
        ones(neco) .* (length(opts.echo_times) > 1 ? parse(Float64, opts.echo_times[2]) : 1.0)
    else
        parse_numbers(opts.echo_times)
    end
    if 1 < length(TEs) == neco
        TEs = TEs[echoes]
    end
    return TEs
end

function get_phase_offset_smoothing_sigma(opts::RomeoOptions)
    if isempty(opts.phase_offset_smoothing_sigma_mm)
        return [7, 7, 7]
    end
    return parse_numbers(opts.phase_offset_smoothing_sigma_mm)
end

function phase_offset_correction!(opts::RomeoOptions, phase, mag, hdr)
    if opts.phase_offset_correction in ("monopolar", "bipolar")
        return run_phase_offset_correction!(opts, phase, mag, hdr)
    end
    size(phase, 5) == 1 || error("multi-channel data requires phase offset correction, --phase-offset-correction off is not possible")
    return as4d(phase), as4d(mag)
end

function run_phase_offset_correction!(opts::RomeoOptions, phase, mag, hdr)
    polarity = opts.phase_offset_correction
    bipolar_correction = polarity == "bipolar"
    neco = opts.number_of_echoes

    TEs = getTEs(opts, neco, collect(1:neco)) # get all TEs here (not only selected)
    if neco != length(TEs) error("Phase offset determination requires all echo times! ($(length(TEs)) given, $neco required)") end

    opts.verbose && CLI.info("Perform phase offset correction with MCPC-3D-S ($polarity)")
    opts.verbose && opts.multi_channel && CLI.info("Perform coil combination with MCPC-3D-S ($polarity)")

    po = zeros(Float32, (size3(phase)..., size(phase, 5)))
    sigma_mm = get_phase_offset_smoothing_sigma(opts)
    sigma_vox = sigma_mm ./ (hdr.pixdim[2], hdr.pixdim[3], hdr.pixdim[4])

    magnitude = mag === nothing ? ones(Float32, size(phase)) : mag
    combined, mcomb = mcpc3ds_typed(phase, magnitude, TEs, po, bipolar_correction, sigma_vox)

    if opts.multi_channel
        opts.verbose && CLI.info("Saving combined_phase, combined_mag and phase_offset")
        save(combined, "combined_phase", opts, hdr)
        save(mcomb, "combined_mag", opts, hdr)
    else
        opts.verbose && CLI.info("Saving corrected_phase and phase_offset")
        save(combined, "corrected_phase", opts, hdr)
    end
    opts.write_phase_offsets && save(po, "phase_offset", opts, hdr)

    mag4 = size(magnitude, 5) != 1 ? mcomb : as4d(mag)
    return combined, mag4
end

# positional arguments are split by type where keyword arguments are not
mcpc3ds_typed(phase, mag, TEs, po, bipolar_correction, sigma) = mcpc3ds(phase, mag; TEs, po, bipolar_correction, sigma)

function get_keyargs(opts::RomeoOptions, mag, mask, TEs, weights)
    opts.verbose && opts.max_seeds != 1 && CLI.info("Maxseeds are $(opts.max_seeds)")
    opts.verbose && opts.merge_regions && CLI.info("Region merging is activated")
    opts.verbose && opts.correct_regions && CLI.info("Region correcting is activated")
    opts.verbose && CLI.info("individual unwrapping is $(opts.individual_unwrapping)")
    opts.verbose && !opts.individual_unwrapping && CLI.info("echo $(opts.template) used as template")
    return (; mag, mask, TEs, weights,
        correctglobal = opts.correct_global,
        maxseeds = opts.max_seeds,
        merge_regions = opts.merge_regions,
        correct_regions = opts.correct_regions,
        wrap_addition = opts.wrap_addition,
        temporal_uncertain_unwrapping = opts.temporal_uncertain_unwrapping,
        individual = opts.individual_unwrapping,
        template = opts.template)
end

function parseweights(opts::RomeoOptions)
    w = opts.weights
    if isfile(w) && any(==('.'), basename(w))
        weights, _ = loadnii(w)
        return UInt8.(as4d(weights))
    end
    flags = ROMEO.parse_weight_flags(w)
    return flags === nothing ? Symbol(w) : flags
end

format_size(a) = "(" * join((string(size(a, i)) for i in 1:ndims(a) if i <= 3 || size(a, i) > 1), ", ") * ")"

function set_mask(opts::RomeoOptions, phase, mag, hdr, TEs, weights)
    m = opts.mask[1]
    if isfile(m)
        opts.verbose && CLI.info("Trying to read mask from file $m")
        mask = first(loadnii(m)) .!= 0
        if size(mask, 4) != 1 || size(mask, 5) != 1 || size3(mask) != size3(phase)
            error("size of mask is $(format_size(mask)), but it should be $(format_size(phase))!")
        end
        return reshape(mask, size3(mask))
    elseif m == "robustmask" && mag !== nothing
        opts.verbose && CLI.info("Calculate robustmask from magnitude, saved as mask.nii")
        template_echo = min(opts.template, size(mag, 4))
        mask = robustmask(mag[:,:,:,template_echo])
        save(mask, "mask", opts, hdr)
        return mask
    elseif contains(m, "qualitymask")
        threshold = if length(opts.mask) > 1
            Float64(parse(Float32, opts.mask[2]))
        elseif length(split(m)) > 1
            Float64(parse(Float32, split(m)[2]))
        else
            0.1 # default threshold
        end
        qmap = voxelquality(phase; get_keyargs(opts, mag, nothing, TEs, weights)...)
        mask = robustmask(qmap; threshold)
        save(mask, "mask", opts, hdr)
        return mask
    elseif m != "nomask"
        error("masking option '$m' is undefined" * (tryparse(Float32, m) isa Float32 ? " (Maybe '-k qualitymask $m' was meant?)" : ""))
    end
    return nothing
end

# The magnitude, the mask and the weights each have several possible types.
# They are split one at a time, so that every call below them has static types
# and a compiled program does not need dynamic dispatch.
function unwrap_and_write(opts::RomeoOptions, phase, mag, hdr, TEs)
    weights = parseweights(opts)
    if mag === nothing
        unwrap_and_write(opts, phase, nothing, weights, hdr, TEs)
    else
        unwrap_and_write(opts, phase, mag, weights, hdr, TEs)
    end
end

function unwrap_and_write(opts::RomeoOptions, phase, mag, weights, hdr, TEs)
    mask = set_mask(opts, phase, mag, hdr, TEs, weights)
    if mask === nothing
        run_unwrap(opts, phase, mag, nothing, weights, hdr, TEs)
    else
        run_unwrap(opts, phase, mag, mask, weights, hdr, TEs)
    end
end

function run_unwrap(opts::RomeoOptions, phase, mag, mask, weights, hdr, TEs)
    keyargs = get_keyargs(opts, mag, mask, TEs, weights)

    opts.verbose && CLI.info("perform unwrapping...")
    regions = zeros(UInt8, size3(phase)) # regions is an output
    unwrap!(phase; keyargs..., regions)
    opts.verbose && CLI.info("unwrapping finished!")

    if opts.max_seeds > 1
        opts.verbose && CLI.info("writing regions...")
        save(regions, "regions", opts, hdr)
    end

    if opts.threshold != Inf
        max = opts.threshold * 2π
        phase[phase .> max] .= 0
        phase[phase .< -max] .= 0
    end

    if opts.mask_unwrapped && mask !== nothing
        phase .*= mask
    end

    save(phase, opts.filename, opts, hdr)

    if !isempty(opts.compute_B0)
        computeB0(opts, phase, mag, hdr, TEs)
    end

    write_qualitymap(opts, phase, keyargs, hdr)
end

function computeB0(opts::RomeoOptions, phase, mag, hdr, TEs)
    if isempty(opts.echo_times)
        error("echo times are required for B0 calculation! Unwrapping has been performed")
    end
    if mag === nothing
        if length(TEs) > 1
            CLI.warn("B0 frequency estimation without magnitude might result in poor handling of noise in later echoes!")
        end
        mag = to_dim(exp.(-TEs/20), Val(4)) # T2*=20ms decay (low value to reduce noise contribution of later echoes)
    end
    B0, snr = compute_B0(phase, mag, TEs, opts.B0_phase_weighting)
    save(B0, opts.compute_B0, opts, hdr)
    save(snr, opts.compute_B0 * "_snr", opts, hdr)
end

# the weighting as a Val, so that every call is statically resolvable
function compute_B0(phase, mag, TEs, weighting::String)
    weighting == "phase_snr" && return B0_and_snr(phase, mag, TEs, Val(:phase_snr))
    weighting == "phase_var" && return B0_and_snr(phase, mag, TEs, Val(:phase_var))
    weighting == "average" && return B0_and_snr(phase, mag, TEs, Val(:average))
    weighting == "TEs" && return B0_and_snr(phase, mag, TEs, Val(:TEs))
    weighting == "mag" && return B0_and_snr(phase, mag, TEs, Val(:mag))
    weighting == "simulated_mag" && return B0_and_snr(phase, mag, TEs, Val(:simulated_mag))
    error("The phase weighting option '$weighting' is not defined!")
end
B0_and_snr(phase, mag, TEs, type::Val) = calculateB0_unwrapped(phase, mag, TEs, type), get_B0_snr(mag, TEs, type)

function write_qualitymap(opts::RomeoOptions, phase, keyargs, hdr)
    # no mask used for writing quality maps
    if opts.write_quality
        opts.verbose && CLI.info("Calculate and write quality map...")
        save(voxelquality(phase; keyargs...), "quality", opts, hdr)
    end
    if opts.write_quality_all
        for i in 1:6
            flags = falses(6)
            flags[i] = true
            opts.verbose && CLI.info("Calculate and write quality map $i...")
            qm = voxelquality(phase; keyargs..., weights=flags)
            if all(qm[1:end-1,1:end-1,1:end-1] .== 1.0)
                opts.verbose && CLI.info("quality map $i skipped for the given inputs")
            else
                save(qm, "quality_$i", opts, hdr)
            end
        end
    end
end

save(image, name, opts::RomeoOptions, hdr) = savenii(image, name, opts.output, hdr)

function saveconfiguration(opts::RomeoOptions, args, version)
    # Cite only what this run used. Phase-offset correction is already resolved
    # by load_data_and_resolve_args! at this point, so "off" here really means the
    # method did not run - including the case where -B switched it on and equal
    # echo times switched it back off.
    cite = [:romeo]
    if opts.multi_channel || opts.phase_offset_correction != "off"
        push!(cite, :mcpc3ds)
    end
    if opts.weights == "bestpath"
        push!(cite, :bestpath)
    end

    write_provenance(opts.output, "romeo";
        version, args, settings = settings_record(opts), cite,
        optional = [:phase_based_masking, :qsmxt, :julia],
        inputs = Pair{String,Union{Nothing,String}}["phase" => opts.phase, "magnitude" => (isempty(opts.magnitude) ? nothing : opts.magnitude)],
        packages = [ROMEO, MriResearchTools],
        describe = describe_input,
    )
end

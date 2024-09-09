function ROMEO.unwrapping_main(args; version="App 4.0")
    settings = getargs(args, version)
    data = load_data_and_resolve_args!(settings)

    mkpath(settings["output"])    
    saveconfiguration(settings["output"], settings, args, version)
    
    ## Perform phase offset correction with MCPC-3D-S (and in case of 5D coil combinadion)
    if settings["phase-offset-correction"] in ["monopolar", "bipolar"]
        phase_offset_correction!(settings, data)
    end
    
    select_echoes!(data, settings)

    set_mask!(data, settings)

    keyargs = get_keyargs(settings, data)

    unwrapping(data, settings, keyargs)
    
    if settings["threshold"] != Inf
        max = settings["threshold"] * 2π
        data["phase"][data["phase"] .> max] .= 0
        data["phase"][data["phase"] .< -max] .= 0
    end

    if settings["mask-unwrapped"] && haskey(data, "mask")
        data["phase"] .*= data["mask"]
    end

    save(data["phase"], settings["filename"], settings)

    if !isempty(settings["compute-B0"])
        computeB0(settings, data)
    end

    write_qualitymap(settings, data, keyargs)

    return 0
end

function load_data_and_resolve_args!(settings)
    settings["filename"] = "unwrapped"
    if endswith(settings["output"], ".nii") || endswith(settings["output"], ".nii.gz")
        settings["filename"] = basename(settings["output"])
        settings["output"] = dirname(settings["output"])
    end

    if settings["weights"] == "romeo"
        if isnothing(settings["magnitude"])
            settings["weights"] = "romeo4"
        else
            settings["weights"] = "romeo3"
        end
    end

    if settings["mask-unwrapped"] && settings["mask"][1] == "nomask"
        settings["mask"][1] = "robustmask"
    end

    if settings["mask"][1] == "robustmask" && isnothing(settings["magnitude"])
        settings["mask"][1] = "nomask"
        @warn "robustmask was chosen but no magnitude is available. No mask is used!" maxlog=2
    end

    settings["mmap-phase"] = !settings["no-mmap"] && !endswith(settings["phase"], ".gz")
    settings["mmap-mag"] = !settings["no-mmap"] && (isnothing(settings["magnitude"]) || !endswith(settings["magnitude"], ".gz"))

    data = Dict{String, AbstractArray}()
    data["phase"] = readphase(settings["phase"], mmap=settings["mmap-phase"], rescale=!settings["no-rescale"])
    settings["verbose"] && println("Phase loaded!")
    if !isnothing(settings["magnitude"])
        data["mag"] = readmag(settings["magnitude"], mmap=settings["mmap-mag"])
        settings["verbose"] && println("Mag loaded!")
    end
    
    settings["header"] = header(data["phase"])
    settings["neco"] = size(data["phase"], 4)

    # activate phase-offset-correction as default (monopolar)
    settings["multi-channel"] = size(data["phase"], 5) > 1
    if (!isempty(settings["compute-B0"]) || settings["multi-channel"] || settings["phase-offset-correction"] == "on") && settings["phase-offset-correction"] ∉ ["bipolar", "off"]
        settings["phase-offset-correction"] = "monopolar"
        settings["verbose"] && println("Phase offset correction with MCPC-3D-S set to monopolar")
    end
    if settings["neco"] == 1
        settings["phase-offset-correction"] = "off"
        settings["verbose"] && println("Phase offset correction with MCPC-3D-S turned off (only one echo)")
    end
    if settings["phase-offset-correction"] == "default off"
        settings["phase-offset-correction"] = "off"
    end

    ## Echoes for unwrapping
    settings["echoes"] = try
        getechoes(settings, settings["neco"])
    catch y
        if isa(y, BoundsError)
            error("echoes=$(join(settings["unwrap-echoes"], " ")): specified echo out of range! Number of echoes is $(settings["neco"])")
        else
            error("echoes=$(join(settings["unwrap-echoes"], " ")) wrongly formatted!")
        end
    end
    settings["verbose"] && println("Echoes are $(settings["echoes"])")

    settings["TEs"] = getTEs(settings, settings["neco"], settings["echoes"])
    settings["verbose"] && println("TEs are $(settings["TEs"])")

    if 1 < length(settings["echoes"]) && length(settings["echoes"]) != length(settings["TEs"])
        error("Number of chosen echoes is $(length(settings["echoes"])) ($(settings["neco"]) in .nii data), but $(length(settings["TEs"])) TEs were specified!")
    end
    
    if haskey(data, "mag") && (size.(Ref(data["mag"]), 1:3) != size.(Ref(data["phase"]), 1:3) || size(data["mag"], 4) < maximum(settings["echoes"]))
        error("size of magnitude and phase does not match!")
    end

    equal_echo_time = length(settings["TEs"]) >= 2 && settings["TEs"][1] == settings["TEs"][2]
    if settings["phase-offset-correction"] != "off" && equal_echo_time
        @warn "The echo times 1 and 2 ($(settings["TEs"])) need to be different for MCPC-3D-S phase offset correction! No phase offset correction performed"
        settings["phase-offset-correction"] = "off"
    end

    return data
end

function phase_offset_correction!(settings, data)
    polarity = settings["phase-offset-correction"]
    bipolar_correction = polarity == "bipolar"

    TEs = getTEs(settings, settings["neco"], :) # get all TEs here (not only selected)
    if settings["neco"] != length(TEs) error("Phase offset determination requires all echo times! ($(length(TEs)) given, $(settings["neco"]) required)") end
    
    settings["verbose"] && println("Perform phase offset correction with MCPC-3D-S ($polarity)")
    settings["verbose"] && settings["multi-channel"] && println("Perform coil combination with MCPC-3D-S ($polarity)")

    po = zeros(eltype(data["phase"]), (size(data["phase"])[1:3]...,size(data["phase"],5)))
    sigma_mm = get_phase_offset_smoothing_sigma(settings)
    sigma_vox = sigma_mm ./ header(data["phase"]).pixdim[2:4]
    
    mag = if haskey(data, "mag") data["mag"] else ones(size(data["phase"])) end
    phase, mcomb = mcpc3ds(data["phase"], mag; TEs, po, bipolar_correction, sigma=sigma_vox)
    data["phase"] = phase
    
    if size(mag, 5) != 1
        data["mag"] = mcomb
    end
    if settings["multi-channel"]
        settings["verbose"] && println("Saving combined_phase, combined_mag and phase_offset")
        save(phase, "combined_phase", settings)
        save(mcomb, "combined_mag", settings)
    else
        settings["verbose"] && println("Saving corrected_phase and phase_offset")
        save(phase, "corrected_phase", settings)
    end
    settings["write-phase-offsets"] && save(po, "phase_offset", settings)
end

function get_keyargs(settings, data)
    keyargs = Dict{Symbol, Any}()

    if haskey(data, "mag")
        keyargs[:mag] = data["mag"]
    end
    if haskey(data, "mask")
        keyargs[:mask] = data["mask"]
    end

    keyargs[:TEs] = settings["TEs"]
    keyargs[:correctglobal] = settings["correct-global"]
    keyargs[:weights] = parseweights(settings)
    keyargs[:maxseeds] = settings["max-seeds"]
    settings["verbose"] && keyargs[:maxseeds] != 1 && println("Maxseeds are $(keyargs[:maxseeds])")
    keyargs[:merge_regions] = settings["merge-regions"]
    settings["verbose"] && keyargs[:merge_regions] && println("Region merging is activated")
    keyargs[:correct_regions] = settings["correct-regions"]
    settings["verbose"] && keyargs[:correct_regions] && println("Region correcting is activated")
    keyargs[:wrap_addition] = settings["wrap-addition"]
    keyargs[:temporal_uncertain_unwrapping] = settings["temporal-uncertain-unwrapping"]
    keyargs[:individual] = settings["individual-unwrapping"]
    settings["verbose"] && println("individual unwrapping is $(keyargs[:individual])")
    keyargs[:template] = settings["template"]
    settings["verbose"] && !settings["individual-unwrapping"] && println("echo $(keyargs[:template]) used as template")

    return keyargs
end

function select_echoes!(data, settings)
    data["phase"] = data["phase"][:,:,:,settings["echoes"]]
    if haskey(data, "mag")
        data["mag"] = data["mag"][:,:,:,settings["echoes"]]
    end
end

function set_mask!(data, settings)
    if isfile(settings["mask"][1])
        settings["verbose"] && println("Trying to read mask from file $(settings["mask"][1])")
        data["mask"] = niread(settings["mask"][1]).raw .!= 0
        if size(data["mask"]) != size(data["phase"])[1:3]
            error("size of mask is $(size(data["mask"])), but it should be $(size(data["phase"])[1:3])!")
        end
    elseif settings["mask"][1] == "robustmask" && haskey(data, "mag")
        settings["verbose"] && println("Calculate robustmask from magnitude, saved as mask.nii")
        template_echo = min(settings["template"], size(data["mag"], 4))
        data["mask"] = robustmask(data["mag"][:,:,:,template_echo])
        save(data["mask"], "mask", settings)
    elseif contains(settings["mask"][1], "qualitymask")
        threshold = if length(settings["mask"]) > 1
            parse(Float32, settings["mask"][2])
        elseif length(split(settings["mask"][1])) > 1
            parse(Float32, split(settings["mask"][1])[2])
        else
            0.1 # default threshold
        end
        qmap = voxelquality(data["phase"]; get_keyargs(settings, data)...)
        data["mask"] = robustmask(qmap; threshold)
        save(data["mask"], "mask", settings)
    elseif settings["mask"][1] != "nomask"
        opt = settings["mask"][1]
        error("masking option '$opt' is undefined" * ifelse(tryparse(Float32, opt) isa Float32, " (Maybe '-k qualitymask $opt' was meant?)", ""))
    end
end

function unwrapping(data, settings, keyargs)
    settings["verbose"] && println("perform unwrapping...")
    regions = zeros(UInt8, size(data["phase"])[1:3]) # regions is an output
    unwrap!(data["phase"]; keyargs..., regions)
    settings["verbose"] && println("unwrapping finished!")

    if settings["max-seeds"] > 1
        settings["verbose"] && println("writing regions...")
        save(regions, "regions", settings)
    end
end

function computeB0(settings, data)
    if isempty(settings["echo-times"])
        error("echo times are required for B0 calculation! Unwrapping has been performed")
    end
    if !haskey(data, "mag")
        if length(settings["TEs"]) > 1
            @warn "B0 frequency estimation without magnitude might result in poor handling of noise in later echoes!"
        end
        data["mag"] = to_dim(exp.(-settings["TEs"]/20), 4) # T2*=20ms decay (low value to reduce noise contribution of later echoes)
    end
    B0 = calculateB0_unwrapped(data["phase"], data["mag"], settings["TEs"])
    save(B0, settings["compute-B0"], settings)
end

function write_qualitymap(settings, data, keyargs)
    # no mask used for writing quality maps
    if settings["write-quality"]
        settings["verbose"] && println("Calculate and write quality map...")
        save(voxelquality(data["phase"]; keyargs...), "quality", settings)
    end
    if settings["write-quality-all"]
        for i in 1:6
            flags = falses(6)
            flags[i] = true
            settings["verbose"] && println("Calculate and write quality map $i...")
            qm = voxelquality(data["phase"]; keyargs..., weights=flags)
            if all(qm[1:end-1,1:end-1,1:end-1] .== 1.0)
                settings["verbose"] && println("quality map $i skipped for the given inputs")
            else
                save(qm, "quality_$i", settings)
            end
        end
    end
end

save(image, name, settings::Dict) = savenii(image, name, settings["output"], settings["header"])

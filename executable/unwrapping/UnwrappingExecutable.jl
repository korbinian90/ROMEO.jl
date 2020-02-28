include("../../romeo/romeo.jl")
include("../../MRI/src/NIfTI_mod.jl")
include("../../MRI/src/utility.jl")
include("../../MRI/src/smoothing.jl")

using ArgParse

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    ret = unwrapping_main(ARGS)
    if ret == 0
        return 0
    else
        println(ret)
        return 1
    end
end

function getargs(args)
    if isempty(args) args = ["--help"] end
    s = ArgParseSettings()
    @add_arg_table s begin
        "phase"
            help = "The phase image used for unwrapping"
        "--magnitude", "-m"
            help = "The magnitude image (better unwrapping if specified)"
        "--output", "-o"
            help = "The output path and filename"
            default = "unwrapped.nii"
        "--echo-times", "-t"
            help = """The relative echo times required for temporal unwrapping (default is 1:n)
                    specified in array or range syntax (eg. [1.5,3.0] or 2:5 or "ones(nr_of_time_points)")"""
        "--mask", "-k"
            help = "<mask_file> | nomask | robustmask"
            default = "robustmask"
        "--individual-unwrapping", "-i"
            help = """Unwraps the echoes individually (not temporal)
                    Temporal unwrapping only works if phase evolves approximatively linear with time"""
            action = :store_true
        "--unwrap-echoes", "-e"
            help = "Unwrap only the specified echoes"
            default = ":"
        "--weights", "-w"
            help = """<4d-weights-file> | romeo3 | romeo2 | bestpath
            romeo2 - magnitude and phase spatial coherence weights are calculated, default for single echo data
            romeo3 - same as romeo2 plus phase temporal coherence weight, default for multi echo data"""
            default = "romeo3"
        "--compute-B0", "-B"
            help = "Calculate combined B0 map in [rad/s]"
            action = :store_true
        "--template", "-T"
            help = "template volume number for template unwrapping, default is 2"

    end
    parse_args(args, s)
end

function unwrapping_main(args)
    settings = getargs(args)

    writedir = settings["output"]
    filename = "unwrapped"
    if occursin(r"\.nii$", writedir)
        filename = basename(writedir)
        writedir = dirname(writedir)
    end

    mkpath(writedir)

    phasenii = readphase(settings["phase"])
    neco = size(phasenii, 4)

    echoes = try
        getechoes(settings, neco)
    catch y
        if isa(y, BoundsError)
            return "ERROR: echoes=$(settings["unwrap-echoes"]): specified echo out of range! Number of echoes is $neco"
        else
            return "ERROR: echoes=$(settings["unwrap-echoes"]) wrongly formatted!"
        end
    end

    hdr = copy(phasenii.header)
    hdr.scl_slope = 1
    hdr.scl_inter = 0

    phase = createniiforwriting(view(phasenii,:,:,:,echoes), filename, writedir; header = hdr, datatype = Float32)

    keyargs = Dict()
    if settings["magnitude"] != nothing
        keyargs[:mag] = view(readmag(settings["magnitude"]).raw,:,:,:,echoes)
        if size(keyargs[:mag]) != size(phase)
            return "ERROR: size of magnitude and phase does not match!"
        end
    end

    ## get settings
    if isfile(settings["mask"])
        keyargs[:mask] = niread(settings["mask"]) .!= 0
        if size(keyargs[:mask]) != size(phase)[1:3]
            return "ERROR: size of mask is $(size(keyargs[:mask])), but it should be $(size(phase)[1:3])!"
        end
    elseif settings["mask"] == "robustmask" && haskey(keyargs, :mag)
        keyargs[:mask] = getrobustmask(keyargs[:mag][:,:,:,1])
        savenii(keyargs[:mask], "mask", writedir; header = hdr)
    end
    if length(echoes) > 1
        keyargs[:TEs] = getTEs(settings, neco, echoes)
    end
    if isfile(settings["weights"]) && splitext(settings["weights"])[2] != ""
        keyargs[:weights] = UInt8.(niread(settings["weights"]))
    else
        keyargs[:weights] = Symbol(settings["weights"])
    end

    if settings["template"] != nothing
        keyargs[:template] = eval(Meta.parse(settings["template"]))
    else
        keyargs[:template] = 2
    end

    ## Error messages
    if 1 < length(echoes) && length(echoes) != length(keyargs[:TEs])
        return "ERROR: Number of chosen echoes is $(length(echoes)) ($neco in .nii data), but $(length(keyargs[:TEs])) TEs were specified!"
    end

    if settings["individual-unwrapping"] && length(echoes) > 1
        kunwrap_single!(phase; keyargs...)
    else
        kunwrap!(phase; keyargs...)
    end

    if settings["compute-B0"]
        if settings["echo-times"] == nothing
            return "ERROR: echo times are required for B0 calculation!"
        end
        if !haskey(keyargs, :mag)
            keyargs[:mag] = ones(1,1,1,size(phase,4))
        end
        TEs = reshape(keyargs[:TEs],1,1,1,:)
        B0 = 1000 * sum(phase ./ TEs .* keyargs[:mag]; dims = 4)
        B0 ./= sum(keyargs[:mag]; dims = 4)

        savenii(B0, "B0", writedir; header = hdr)
    end

    @show writedir
    return 0
end

function getechoes(settings, neco)
    echoes = eval(Meta.parse(settings["unwrap-echoes"]))
    if typeof(echoes) <: Int
        echoes = [echoes]
    end
    echoes = (1:neco)[echoes]
    if length(echoes) == 1 echoes = echoes[1] end
    return echoes
end

function getTEs(settings, neco, echoes)
    if settings["echo-times"] != nothing
        TEs = eval(Meta.parse(settings["echo-times"]))
        if length(TEs) == neco
            TEs = TEs[echoes]
        end
    else
        TEs = (1:neco)[echoes]
    end
    return TEs
end

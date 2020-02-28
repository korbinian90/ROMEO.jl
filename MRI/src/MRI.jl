module MRI

include("NIfTI_mod.jl")
include("utility.jl")
include("smoothing.jl")
include("laplacianunwrapping.jl")
include("intensitycorrection.jl")
include("VSMbasedunwarping.jl")

export Data,
        readphase, readmag, niread,
        savenii,
        createniiforwriting,
        getrobustmask, robustmask!,
        robustrescale,
        combine_echoes,
        getHIP,
        laplacianunwrap, laplacianunwrap!,
        getVSM,
        unwarp,
        thresholdforward,
        gaussiansmooth3d!, gaussiansmooth3d,
        makehomogeneous!, makehomogeneous,
        getsensitivity,
        getscaledimage,
        estimatequantile,
        RSS

end # module

# TODO create tests

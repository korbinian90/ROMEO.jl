function unwrap!(wrapped::AbstractArray{T, 3}; weights = :romeo, keyargs...) where {T <: AbstractFloat}
    nbins = 256

    weights = calculateweights(wrapped, weights, nbins; keyargs...)
    @assert sum(weights) != 0 "Unwrap-weights are all zero!"

    seed = findseed(weights)
    if haskey(keyargs, :phase2) && haskey(keyargs, :TEs) # requires multiecho
        seedcorrection!(wrapped, seed, keyargs[:phase2], keyargs[:TEs])
    end

    growRegionUnwrap!(wrapped, weights, seed, nbins, nothing)
end

function unwrap4d!(wrapped::AbstractArray{T, 4}; weights=:romeo, keyargs...) where {T <: AbstractFloat}
    nbins = 256
    ref = size(wrapped, 4) >= 3 ? 3 : 1
    template = 2
    TEs = keyargs[:TEs]
    @show typeof(keyargs[:mag]) typeof(TEs)
    data = Data(keyargs[:mag], zeros(Float32, size(wrapped)[1:3]), zeros(Float32, size(wrapped)[1:3]), TEs)

    args = Dict{Symbol, Any}(keyargs)
    args[:phase2] = wrapped[:,:,:,ref]
    args[:TEs] = TEs[[template, ref]]
    if haskey(args, :mag)
        args[:mag] = args[:mag][:,:,:,template]
    end
    weights = calculateweights(wrapped[:,:,:,template], weights, nbins; args...)
    @assert sum(weights) != 0 "Unwrap-weights are all zero!"

    seed = findseed(weights)
    if haskey(keyargs, :phase2) && haskey(keyargs, :TEs) # requires multiecho
        seedcorrection!(wrapped, seed, keyargs[:phase2], keyargs[:TEs])
    end

    growRegionUnwrap!(wrapped, weights, seed, nbins, data), data
end

struct Data
    mag::Array{Float32,4}
    B0::Array{Float32,3}
    PO::Array{Float32,3}
    TEs::Array{Float32,1}
end

"""
    unwrap(wrapped::AbstractArray{T, 3}; weights = :romeo, keyargs...)

ROMEO unwrapping.
Options for weights are :romeo, :bestpath.

###  Optional keyword arguments:

- `mag`: Additional mag weights are used.
- `mask`: Unwrapping is only performed inside the mask.
- `phase2`: A second reference phase image (possibly with different echo time).
   It is used for calculating the phasecoherence weight.
- `TEs`: The echo times of the phase and the phase2 images.

"""
unwrap(wrapped; keyargs...) = unwrap!(copy(wrapped); keyargs...)

unwrap!(wrapped::AbstractArray{T,2}; keyargs...) where {T <: AbstractFloat} = unwrap!(reshape(wrapped, size(wrapped)..., 1); keyargs...)
unwrap!(wrapped::AbstractArray{T,1}; keyargs...) where {T <: AbstractFloat} = unwrap!(reshape(wrapped, size(wrapped)..., 1); keyargs...)

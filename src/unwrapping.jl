function unwrap!(wrapped::AbstractArray{T, 3}; weights = :romeo, keyargs...) where {T <: AbstractFloat}
    nbins = 256

    weights = calculateweights(wrapped, weights, nbins; keyargs...)
    @assert sum(weights) != 0 "Unwrap-weights are all zero!"

    seed = findseed(wrapped, weights)
    if haskey(keyargs, :phase2) && haskey(keyargs, :TEs) # requires multiecho
        seedcorrection!(wrapped, seed, keyargs[:phase2], keyargs[:TEs])
    end

    growRegionUnwrap!(wrapped, weights, seed, nbins)
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

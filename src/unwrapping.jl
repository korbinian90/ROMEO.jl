function unwrap!(wrapped::AbstractArray{T, 3}; weights = :romeo, keyargs...) where {T <: AbstractFloat}
    nbins = 256

    weights = calculateweights(wrapped, weights, nbins; keyargs...)
    @assert sum(weights) != 0 "Unwrap-weights are all zero!"

    growRegionUnwrap!(wrapped, weights, nbins; keyargs...)

    if haskey(keyargs, :correctglobal) && keyargs[:correctglobal]
        mask = if haskey(keyargs, :mask)
            keyargs[:mask]
        else
            trues(size(wrapped))
        end
        wrapped .-= (2π * median(round.(wrapped[mask] ./ 2π)))
    end
    return wrapped
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
- correctglobal: if true corrects for global n2π offsets

"""
unwrap(wrapped; keyargs...) = unwrap!(copy(wrapped); keyargs...)

unwrap!(wrapped::AbstractArray{T,2}; keyargs...) where {T <: AbstractFloat} = unwrap!(reshape(wrapped, size(wrapped)..., 1); keyargs...)
unwrap!(wrapped::AbstractArray{T,1}; keyargs...) where {T <: AbstractFloat} = unwrap!(reshape(wrapped, size(wrapped)..., 1); keyargs...)

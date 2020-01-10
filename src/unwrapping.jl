function unwrap!(wrapped::AbstractArray{T, 3}; weights = :romeo, keyargs...) where {T <: AbstractFloat}
    nbins = 256

    weights = calculateweights(wrapped, nbins, weights; keyargs...)
    @assert sum(weights) != 0 "Unwrap-weights are all zero!"

    seed = findseed(wrapped, weights)
    if haskey(keyargs, :phase2) && haskey(keyargs, :TEs) # requires multiecho
        seedcorrection!(wrapped, seed, keyargs[:phase2], keyargs[:TEs])
    end

    growRegionUnwrap!(wrapped, weights, seed, nbins)
end

# unwrap version that does not modify its input
unwrap(wrapped; keyargs...) = unwrap!(copy(wrapped); keyargs...)

unwrap!(wrapped::AbstractArray{T,2}; keyargs...) where {T <: AbstractFloat} = unwrap!(reshape(wrapped, size(wrapped)..., 1); keyargs...)
unwrap!(wrapped::AbstractArray{T,1}; keyargs...) where {T <: AbstractFloat} = unwrap!(reshape(wrapped, size(wrapped)..., 1); keyargs...)

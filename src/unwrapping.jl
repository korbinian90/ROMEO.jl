function unwrap!(wrapped::AbstractArray{T,3}; weights=:romeo, keyargs...) where {T <: AbstractFloat}
    nbins = 256

    weights = calculateweights(wrapped, weights, nbins; keyargs...)
    @assert sum(weights) != 0 "Unwrap-weights are all zero!"

    seed = findseed(wrapped, weights)
    if haskey(keyargs, :phase2) && haskey(keyargs, :TEs) # requires multiecho
        seedcorrection!(wrapped, seed, keyargs[:phase2], keyargs[:TEs])
    end

    growRegionUnwrap!(wrapped, weights, seed, nbins)

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

function unwrap!(
    wrapped::AbstractArray{T, 4}; TEs=1:size(wrapped, 4), template=2, p2ref=1, keyargs...
) where {T <: AbstractFloat}
    args = Dict{Symbol, Any}(keyargs)
    args[:phase2] = wrapped[:,:,:,p2ref]
    args[:TEs] = TEs[[template, p2ref]]
    if haskey(args, :mag)
        args[:mag] = args[:mag][:,:,:,template]
    end
    unwrap!(view(wrapped,:,:,:,template); args...)
    for ieco in [(template-1):-1:1; (template+1):length(TEs)]
        iref = if (ieco < template) ieco+1 else ieco-1 end
        refvalue = wrapped[:,:,:,iref] .* (TEs[ieco] / TEs[iref])
        wrapped[:,:,:,ieco] .= unwrapvoxel.(wrapped[:,:,:,ieco], refvalue)
    end
    return wrapped
end

unwrap_individual(wrapped; keyargs...) = unwrap_individual!(Float32.(wrapped); keyargs...)
function unwrap_individual!(
    wrapped::AbstractArray{T, 4}; TEs=1:size(wrapped, 4), keyargs...
) where {T <: AbstractFloat}
    for i in 1:length(TEs)
        e2 = if (i == 1) 2 else i-1 end
        echoes = [i, e2]
        args = Dict()
        if haskey(keyargs, :mag) args[:mag] = keyargs[:mag][:,:,:,i] end
        unwrap!(view(wrapped,:,:,:,i); phase2=wrapped[:,:,:,e2], TEs=TEs[echoes], args...)
    end
    return wrapped
end

function calculateweights(wrapped, nbins, weights; keyargs...)
    if weights == :romeo
        return calculateweights(wrapped, nbins; keyargs...)
    elseif weights == :bestpath
        return calculateweights_bestpath(wrapped, nbins; keyargs...)
    end
end

# each edge gets a weight
function calculateweights(wrapped, nbins; keyargs...)
    if haskey(keyargs, :mag)
        args = Dict{Symbol, Any}(keyargs)
        args[:maxmag] = maximum(args[:mag])
        if haskey(keyargs, :mask)
            args[:mag] = Float32.(args[:mag]) .* args[:mask]
        end
        keyargs = NamedTuple{Tuple(keys(args))}(values(args))
    end
    mask = if haskey(keyargs, :mask)
        keyargs[:mask]
    else
        trues(size(wrapped))
    end

    stridelist = strides(wrapped)
    weights = zeros(UInt8, 3, size(wrapped)...)

    for dim in 1:3
        neighbor = stridelist[dim]
        for I in LinearIndices(wrapped)
            if mask[I]
                weights[dim + (I-1)*3] = getweight(wrapped, I, neighbor, nbins, dim; keyargs...)
            end
        end
    end
    return weights
end

## best path weights
# Abdul-Rahamn https://doi.org/10.1364/AO.46.006623

function calculateweights_bestpath(wrapped, nbins; keyargs...)
    scale(w) = UInt8.(min(max(round((1 - (w / 10)) * (nbins - 1)), 1), 255)) # scaling function
    weights = scale.(getbestpathweight(wrapped))
    if haskey(keyargs, :mask) # apply mask to weights
        mask = keyargs[:mask]
        weights .*= reshape(mask, 1, size(mask)...)
    end
end

function getbestpathweight(φ)
    R = 1 ./ getD(φ)
    weight = zeros(3, size(R)...)
    for idim in 1:3
        n = strides(R)[idim]
        for i in 1:length(R)-n
            weight[idim + 3i] = R[i] + R[i+n]
        end
    end
    return weight
end

function getD(φ)
    directions = Iterators.product(-1:1,-1:1,-1:1)
    neighbors = unique(abs(sum(strides(φ) .* d)) for d in directions if d != (0,0,0))
    D2 = zeros(size(φ))
    @inbounds for n in neighbors, i in 1+n:length(φ)-n
        D2[i] += (γ(φ[i-n] - φ[i]) - γ(φ[i] - φ[i+n])) ^ 2
    end
    return sqrt.(D2)
end

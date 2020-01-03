# each edge gets a weight
function calculateweights(wrapped, nbins, weights; keyargs...)
    if weights == :romeo
        return calculateweights(wrapped, nbins; keyargs...)
    elseif weights == :bestpath
        return calculateweights_bestpath(wrapped, nbins; keyargs...)
    end
end

## ROMEO weights

function calculateweights(wrapped, nbins; keyargs...)
    mask, keyargs = parsekeyargs(keyargs, wrapped)
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

function parsekeyargs(keyargs, wrapped)
    if haskey(keyargs, :mag)
        args = Dict{Symbol, Any}(keyargs)
        if haskey(keyargs, :mask)
            @show "here"
            args[:mag] = Float32.(args[:mag]) .* args[:mask]
        end
        args[:maxmag] = maximum(args[:mag])
        keyargs = NamedTuple{Tuple(keys(args))}(values(args))
    end
    mask = if haskey(keyargs, :mask)
        keyargs[:mask]
    else
        trues(size(wrapped))
    end
    return mask, keyargs
end

@inline function getweight(P, i, k, nbins, dim; keyargs...) # Phase, index, neighbor-offset, nbins, dim
    j = i+k
    if !checkbounds(Bool, P, j) return 0 end

    phasecoherence = 1 - abs(γ(P[i] - P[j]) / π)

    phasegradientcoherence = 1
    if haskey(keyargs, :phase2)
        P2, TEs = keyargs[:phase2], keyargs[:TEs]
        phasegradientcoherence = max(0, 1 - abs(γ(P[i] - P[j]) - γ(P2[i] - P2[j]) * TEs[1] / TEs[2]))
    end

    weight = phasecoherence * phasegradientcoherence

    if haskey(keyargs, :mag)
        M, globalmaxmag = keyargs[:mag], keyargs[:maxmag]
        mini, maxi = minmax(M[i], M[j])
        magcoherence = (mini / maxi) ^ 2

        magweight = 0.5 + 0.5min(1, mini / (0.5 * globalmaxmag))
        magweight2 = 0.5 + 0.5min(1, (0.5 * globalmaxmag) / maxi) # too high magnitude is not good either (flow artifact)

        weight *= magcoherence * magweight * magweight2
    end

    # weight: 1 is best and 0 worst
    # resaled: 1 is best, nbins is worst, 0 is not valid (not added to queue)
    rescale(w) = max(round(Int, (1 - w) * (nbins - 1)), 1)

    if 0 ≤ weight ≤ 1
        return rescale(weight)
    else
        return 0
    end
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

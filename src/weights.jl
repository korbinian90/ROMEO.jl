# each edge gets a weight
function calculateweights(wrapped, nbins=256, weights=:romeo; keyargs...)
    if weights == :romeo
        return calculateweights_romeo(wrapped, nbins; keyargs...)
    elseif weights == :bestpath
        return calculateweights_bestpath(wrapped, nbins; keyargs...)
    end
end

function calculateweights_romeo(wrapped, nbins; keyargs...)
    mask, P2, TEs, M, maxmag = parsekeyargs(keyargs, wrapped)
    stridelist = strides(wrapped)
    weights = zeros(UInt8, 3, size(wrapped)...)
    # rescale
    ## from: 1 is best and 0 worst
    ## to: 1 is best, nbins is worst, 0 is not valid (not added to queue)
    rescale(w) = if 0 ≤ w ≤ 1
            max(round(Int, (1 - w) * (nbins - 1)), 1)
        else 0
        end

    for dim in 1:3
        neighbor = stridelist[dim]
        for I in LinearIndices(wrapped)
            J = I + neighbor
            if mask[I] && checkbounds(Bool, wrapped, J)
                w = getweight(wrapped, I, J, P2, TEs, M, maxmag)
                weights[dim + (I-1)*3] = rescale(w)
            end
        end
    end
    return weights
end

function parsekeyargs(keyargs, wrapped)
    getval(key) = if haskey(keyargs, key) keyargs[key] else nothing end

    mag = getval(:mag)
    mask = getval(:mask)

    if mask != nothing
        if mag != nothing
            mag *= mask
        end
    else
        mask = trues(size(wrapped))
    end

    maxmag = if mag != nothing
        maximum(mag) else nothing
    end

    return mask, getval(:phase2), getval(:TEs), mag, maxmag
end

phasecoherence(P, i, j) = 1 - abs(γ(P[i] - P[j]) / π)
phasegradientcoherence(P, P2, TEs, i, j) = max(0, 1 - abs(γ(P[i] - P[j]) - γ(P2[i] - P2[j]) * TEs[1] / TEs[2]))
magcoherence(small, big) = (small / big) ^ 2
magweight(small, max) = 0.5 + 0.5min(1, small / (0.5 * max))
magweight2(big, max) = 0.5 + 0.5min(1, (0.5 * max) / big) # too high magnitude is not good either (flow artifact)

# calculates weight of one edge
function getweight(P, i, j, P2, TEs, M, maxmag) # Phase, index, neighbor
    weight = phasecoherence(P, i, j)

    if P2 != nothing && TEs != nothing
        weight *= phasegradientcoherence(P, P2, TEs, i, j)
    end

    if M != nothing && maxmag != nothing
        small, big = minmax(M[i], M[j])
        weight *= magcoherence(small,big) * magweight(small,max) * magweight2(big,max)
    end

    return weight
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

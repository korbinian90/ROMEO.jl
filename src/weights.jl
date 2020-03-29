"""
    calculateweights(wrapped, weights=:romeo, nbins=256; kwargs...)

Calculates weights for all edges.
Options for weights are :romeo, :bestpath.

###  Optional keyword arguments:

- `mag`: Additional mag weights are used.
- `mask`: Unwrapping is only performed inside the mask.
- `phase2`: A second reference phase image (possibly with different echo time).
   It is used for calculating the phasecoherence weight.
- `TEs`: The echo times of the phase and the phase2 images.

"""
function calculateweights(wrapped, weights=:romeo, nbins=256; kwargs...)
    weights = if weights == :romeo
        calculateweights_romeo(wrapped, nbins; kwargs...)
    elseif weights == :bestpath
        calculateweights_bestpath(wrapped, nbins; kwargs...)
    else
        throw(ArgumentError("Weight $weight not defined!"))
    end
    # these edges do not exist
    weights[1,end,:,:] .= 0
    weights[2,:,end,:] .= 0
    weights[3,:,:,end] .= 0
    return weights
end

# rescale
## from: 1 is best and 0 worst
## to: 1 is best, nbins is worst, 0 is not valid (not added to queue)
function rescale(nbins, w)
    if 0 ≤ w ≤ 1
        max(round(Int, (1 - w) * (nbins - 1)), 1)
    else
        0
    end
end

function calculateweights_romeo(wrapped, nbins, ::Type{T}=UInt8; kwargs...) where T
    mask, P2, TEs, M, maxmag = parsekwargs(kwargs, wrapped)
    stridelist = strides(wrapped)
    weights = zeros(T, 3, size(wrapped)...)
    for dim in 1:3
        neighbor = stridelist[dim]
        for I in LinearIndices(wrapped)
            J = I + neighbor
            if mask[I] && checkbounds(Bool, wrapped, J)
                w = getpweight(wrapped, I, J, P2, TEs)
                if M != nothing w *= getmweight(M, I, J, maxmag) end
                weights[dim + (I-1)*3] = rescale(nbins, w)
            end
        end
    end
    return weights
end

# calculates weight of one edge
function getpweight(P, i, j, P2, TEs) # Phase, index, neighbor
    weight = phasecoherence(P, i, j)
    if P2 != nothing && TEs != nothing
        weight * phasegradientcoherence(P, P2, TEs, i, j)
    else
        weight * phaselinearity(P, i, j)
    end
end

function getmweight(M, i, j, maxmag)
    small, big = minmax(M[i], M[j])
    magcoherence(small,big) * magweight(small,maxmag) * magweight2(big,maxmag)
end

phasecoherence(P, i, j) = 1 - abs(γ(P[i] - P[j]) / π)
phasegradientcoherence(P, P2, TEs, i, j) = max(0, 1 - abs(γ(P[i] - P[j]) - γ(P2[i] - P2[j]) * TEs[1] / TEs[2]))
magcoherence(small, big) = (small / big) ^ 2
magweight(small, maxmag) = 0.5 + 0.5min(1, small / (0.5 * maxmag))
magweight2(big, maxmag) = 0.5 + 0.5min(1, (0.5 * maxmag) / big) # too high magnitude is not good either (flow artifact)

phaselinearity(P, i, j, k) = max(0, 1 - abs(rem2pi(P[i] - 2P[j] + P[k], RoundNearest)))
function phaselinearity(P, i, j)
    weight = 1
    h = 2i-j
    if checkbounds(Bool, P, h)
        weight *= phaselinearity(P, h, i, j)
    end
    k = 2j-i
    if checkbounds(Bool, P, k)
        weight *= phaselinearity(P, i, j, k)
    end
    weight
end


function parsekwargs(kwargs, wrapped)
    getval(key) = if haskey(kwargs, key) kwargs[key] else nothing end
    mag = getval(:mag)
    mask = getval(:mask)
    if mask != nothing
        if mag != nothing
            mag .*= mask
        end
    else
        mask = trues(size(wrapped))
    end
    maxmag = if mag != nothing
        maximum(mag[isfinite.(mag)]) else nothing
    end
    return mask, getval(:phase2), getval(:TEs), mag, maxmag
end

## best path weights
# Abdul-Rahamn https://doi.org/10.1364/AO.46.006623

function calculateweights_bestpath(wrapped, nbins; mask=nothing)
    scale(w) = UInt8.(min(max(round((1 - (w / 10)) * (nbins - 1)), 1), 255))
    weights = scale.(getbestpathweight(wrapped))
    if !(mask isa Nothing) # apply mask to weights
        weights .*= reshape(mask, 1, size(mask)...)
    end
    weights
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

## weights
function getweight(P, i, j, P2, TEs, M, maxmag, flags) # Phase, index, neighbor, ...
    weight = 1.0
    if flags[1] weight *= (0.1 + 0.9phasecoherence(P, i, j)) end
    if flags[2] weight *= (0.1 + 0.9phasegradientcoherence(P, P2, TEs, i, j)) end
    if flags[3] weight *= (0.1 + 0.9phaselinearity(P, i, j)) end
    if !isnothing(M)
        small, big = minmax(M[i], M[j])
        if flags[4] weight *= (0.1 + 0.9magcoherence(small, big)) end
        if flags[5] weight *= (0.1 + 0.9magweight(small, maxmag)) end
        if flags[6] weight *= (0.1 + 0.9magweight2(big, maxmag)) end
    end
    return weight
end

phasecoherence(P, i, j) = 1 - abs(γ(P[i] - P[j]) / π)
phasegradientcoherence(P, P2, TEs, i, j) = max(0, 1 - abs(γ(P[i] - P[j]) - γ(P2[i] - P2[j]) * TEs[1] / TEs[2]))
magcoherence(small, big) = (small / big) ^ 2
magweight(small, maxmag) = 0.5 + 0.5min(1, small / (0.5 * maxmag))
magweight2(big, maxmag) = 0.5 + 0.5min(1, (0.5 * maxmag) / big) # too high magnitude is not good either (flow artifact)

function phaselinearity(P, i, j, k)
    pl = max(0, 1 - abs(rem2pi(P[i] - 2P[j] + P[k], RoundNearest)/2))
    if isnan(pl)
        return 0.5
    end
    return pl
end
function phaselinearity(P, i, j)
    neighbor = j - i
    h = i - neighbor
    k = j + neighbor
    if 0 < h && k <= length(P)
        return phaselinearity(P, h, i, j) * phaselinearity(P, i, j, k)
    else
        return 0.9
    end
end

"""
    calculateweights(wrapped; weights=:romeo, kwargs...)

Calculates weights for all edges.
size(weights) == [3, size(wrapped)...]

###  Optional keyword arguments:

- `weights`: Options are [`:romeo`] | `:romeo2` | `:romeo3` | `:bestpath`.
- `mag`: Additional mag weights are used.
- `mask`: Unwrapping is only performed inside the mask.
- `phase2`: A second reference phase image (possibly with different echo time).
   It is used for calculating the phasecoherence weight.
- `TEs`: The echo times of the phase and the phase2 images as a tuple (eg. (5, 10) or [5, 10]).

"""
function calculateweights(wrapped; weights=:romeo, kwargs...)
    weights = if weights isa AbstractArray{<:Real,4}
        weights
    elseif weights == :bestpath
        calculateweights_bestpath(wrapped; kwargs...)
    else
        calculateweights_romeo(wrapped, weights; kwargs...)
    end
    # these edges do not exist
    weights[1,end,:,:] .= 0
    weights[2,:,end,:] .= 0
    weights[3,:,:,end] .= 0
    return weights
end

calculateweights_romeo(wrapped, weights::AbstractArray{T,4}; kw...) where T = UInt8.(weights)
# defines default weights
function calculateweights_romeo(wrapped, weights::Symbol; kwargs...)
    flags = falses(6)
    if weights == :romeo
        flags[1:4] .= true
    elseif weights == :romeo2
        flags[[1,4]] .= true # phasecoherence, magcoherence
    elseif weights == :romeo3
        flags[[1,2,4]] .= true # phasecoherence, phasegradientcoherence, magcoherence
    elseif weights == :romeo4
        flags[1:4] .= true
    elseif weights == :romeo6
        flags[1:6] .= true # additional magweight, magweight2
    else
        throw(ArgumentError("Weight '$weights' not defined!"))
    end
    return calculateweights_romeo(wrapped, flags; kwargs...)
end

function calculateweights_romeo(wrapped, flags::AbstractArray{Bool,1}; type::Type{T}=UInt8, rescale=rescale, kwargs...) where T
    mask, P2, TEs, M, maxmag = parsekwargs(kwargs, wrapped)
    flags = updateflags(flags, P2, TEs, M)
    stridelist = getdimoffsets(wrapped)
    weights = zeros(T, 3, size(wrapped)...)
    Threads.@threads for dim in 1:3
        neighbor = stridelist[dim]
        for I in LinearIndices(wrapped)
            J = I + neighbor
            if mask[I] && checkbounds(Bool, wrapped, J)
                w = getweight(wrapped, I, J, P2, TEs, M, maxmag, flags)
                weights[dim + (I-1)*3] = rescale(w)
            end
        end
    end
    return weights
end

## utility functions
function parsekwargs(kwargs, wrapped)
    getval(key) = if haskey(kwargs, key) kwargs[key] else nothing end
    mag = getval(:mag)
    mask = getval(:mask)
    if !isnothing(mask)
        if !isnothing(mag)
            mag = mag .* mask
        end
    else
        mask = trues(size(wrapped))
    end
    maxmag = if !isnothing(mag)
        quantile(mag[isfinite.(mag)], 0.95) else nothing
    end
    return mask, getval(:phase2), getval(:TEs), mag, maxmag
end

function updateflags(flags, P2, TEs, M)
    flags = copy(flags)
    if isnothing(M)
        flags[4:6] .= false
    end
    if isnothing(P2) || isnothing(TEs)
        flags[2] = false
    end
    return flags
end

# from: 1 is best and 0 worst
# to: 1 is best, NBINS is worst, 0 is not valid (not added to queue)
function rescale(w)
    if 0 ≤ w ≤ 1
        max(round(Int, (1 - w) * (NBINS - 1)), 1)
    else
        0
    end
end

## best path weights
# Abdul-Rahamn https://doi.org/10.1364/AO.46.006623

function calculateweights_bestpath(wrapped; kwargs...)
    scale(w) = UInt8.(min(max(round((1 - (w / 10)) * (NBINS - 1)), 1), 255))
    weights = scale.(getbestpathweight(wrapped))
    if haskey(kwargs, :mask) # apply mask to weights
        mask = kwargs[:mask]
        weights .*= reshape(mask, 1, size(mask)...)
    end
    weights
end

function getbestpathweight(φ)
    R = 1 ./ getD(φ)
    weight = zeros(3, size(R)...)
    for idim in 1:3
        n = getdimoffsets(R)[idim]
        @inbounds for i in 1:length(R)-n
            weight[idim + 3i] = R[i] + R[i+n]
        end
    end
    return weight
end

function getD(φ)
    directions = Iterators.product(-1:1,-1:1,-1:1)
    neighbors = unique(abs(sum(getdimoffsets(φ) .* d)) for d in directions if d != (0,0,0))
    D2 = zeros(size(φ))
    @inbounds for n in neighbors, i in 1+n:length(φ)-n
        D2[i] += (γ(φ[i-n] - φ[i]) - γ(φ[i] - φ[i+n])) ^ 2
    end
    return sqrt.(D2)
end

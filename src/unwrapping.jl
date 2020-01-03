function unwrap!(wrapped::AbstractArray{T, 3}; weights = :romeo, keyargs...) where {T <: AbstractFloat}
    nbins = 256

    weights = calculateweights(wrapped, nbins, weights; keyargs...)
    @assert sum(weights) != 0 "Unwrap weights are all zero!"

    seed = findseed(wrapped, weights)
    if haskey(keyargs, :phase2) && haskey(keyargs, :TEs) # requires multiecho
        seedcorrection!(wrapped, seed, keyargs[:phase2], keyargs[:TEs])
    end

    growRegionUnwrap!(wrapped, weights, seed, nbins)
end

function seedcorrection!(wrapped, seed, phase2, TEs)
    vox = getfirstvoxfromedge(seed)
    best = Inf
    offset = 0
    for off1 in -2:2, off2 in -1:1
        diff = abs((wrapped[vox] + 2π*off1) / TEs[1] - (phase2[vox] + 2π*off2) / TEs[2])
        diff += (abs(off1) + abs(off2)) / 100 # small panelty for wraps (if TE1 == 2*TE2 wrong value is chosen otherwise)
        if diff < best
            best = diff
            offset = off1
        end
    end
    wrapped[vox] += 2π * offset
    return offset
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

    if 0 ≤ weight ≤ 1 # weight of 1 is best and 0 worst
        return max(round(Int, (1 - weight) * (nbins - 1)), 1)
    else
        return 0
    end # 1 is best, nbins is worst, 0 is not valid (not added to queue)
end

function findseed(wrapped, weights)
    cp = copy(weights)
    cp[cp .== 0] .= 255
    filtered = dilate(cp, 2:4)
    (_, ind) = findmin(filtered)
    return LinearIndices(weights)[ind]
end

# unwrap version that does not modify its input
unwrap(wrapped; keyargs...) = unwrap!(Float32.(wrapped); keyargs...)

# multi echo unwrapping
function unwrap!(wrapped::AbstractArray{T, 4}; TEs = 1:size(wrapped, 4), template = 2, p2ref = 1, keyargs...) where {T <: AbstractFloat}
    args = Dict{Symbol, Any}(keyargs)
    args[:phase2] = view(wrapped,:,:,:,p2ref)
    args[:TEs] = TEs[[template, p2ref]]
    if haskey(args, :mag)
        args[:mag] = view(args[:mag],:,:,:,template)
    end
    unwrap!(view(wrapped,:,:,:,template); args...)
    for iEco in [(template-1):-1:1; (template+1):length(TEs)]
        iRef = if (iEco < template) iEco + 1 else iEco - 1 end
        wrapped[:,:,:,iEco] .= unwrapvoxel.(wrapped[:,:,:,iEco], wrapped[:,:,:,iRef] .* (TEs[iEco] / TEs[iRef]))
        #magf = if !haskey(keyargs, :mag) nothing else view(mag,:,:,:,iEco) end
        #unwrapfilter!(view(wrapped,:,:,:,iEco), magf)
    end
    return wrapped
end

function unwrapfilter!(phase, mag)
    smoothedphase = if mag == nothing
        gaussiansmooth3d(phase, boxsizes = [3,3,3], nbox = 1)
    else
        gaussiansmooth3d(phase; weight = Float32.(mag), boxsizes = [3,3,3], nbox = 1) # corresponds to weighted meanfilter with nbox = 1# TODO try other filter
    end
    phase .= unwrapvoxel.(phase, smoothedphase)
end

unwrap_single(wrapped; keyargs...) = unwrap_single!(Float32.(wrapped); keyargs...)
function unwrap_single!(wrapped::AbstractArray{T1, 4}; TEs = 1:size(wrapped, 4), keyargs...) where {T1 <: AbstractFloat}
    for ieco in 1:length(TEs)
        e2 = ieco - 1
        if e2 == 0 e2 = 2 end
        args = Dict()
        if haskey(keyargs, :mag) args[:mag] = view(keyargs[:mag],:,:,:,ieco) end
        unwrap!(view(wrapped,:,:,:,ieco); phase2 = view(wrapped,:,:,:,e2), TEs = TEs[[ieco, e2]], args...)
    end
    return wrapped
end

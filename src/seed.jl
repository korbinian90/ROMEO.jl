function findseed(weights)
    cp = copy(weights)
    cp[cp .== 0] .= 255
    filtered = dilate(cp, 2:4)
    (_, ind) = findmin(filtered)
    return LinearIndices(weights)[ind]
end

function getseedqueue(weights, nbins)
    #cp = copy(weights)
    #cp[cp .== 0] .= 255
    #filtered = dilate(cp, 2:4)
    queue = PQueue{Int}(nbins)
    for (i, w) in enumerate(maximum(weights; dims=1))
        if w > 0
            push!(queue, i, w)
        end
    end
    return queue
end

function findseed!(queue::PQueue, weights, visited)
    while !isempty(queue)
        ind = pop!(queue)
        if visited[ind] == 0
            return ind
        end
    end
    return 0
end

function findminedge(weights, ind)
    edges = 3(ind-1) .+ (1:3) # get 3 edges for this one index
    minedge = edges[1]
    wmin = weights[minedge]
    for e in edges[2:end]
        if 0 < weights[e] < wmin
            minedge = e
            wmin = weights[e]
        end
    end
    return minedge
end

function findseed(weights, visited)
    cp = copy(weights)
    cp[(reshape(visited, 1, size(visited)...) .!= 0) .| (cp .== 0)] .= 255
    (_, ind) = findmin(cp)
    return LinearIndices(weights)[ind]
end

function seedcorrection!(wrapped, vox, keyargs)
    if haskey(keyargs, :phase2) && haskey(keyargs, :TEs) # requires multiecho
        phase2 = keyargs[:phase2]
        TEs = keyargs[:TEs]
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
    else
        wrapped[vox] = rem2pi(wrapped[vox], RoundNearest)
    end
end

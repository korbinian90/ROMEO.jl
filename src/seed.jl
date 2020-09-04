# returns a function that can repeatedly and efficiently create new seeds
function getseedfunction(seeds, pqueue, visited, weights, wrapped, keyargs)
    seedqueue = getseedqueue(weights)
    notvisited(i) = checkbounds(Bool, visited, i) && (visited[i] == 0)
    stridelist = getdimoffsets(wrapped)
    function addseed!()
        seed = findseed!(seedqueue, weights, visited)
        if seed == 0
            return 255
        end
        for i in 1:6 # 6 directions
            e = getnewedge(seed, notvisited, stridelist, i)
            if e != 0 && weights[e] > 0
                push!(pqueue, e, weights[e])
            end
        end
        seedcorrection!(wrapped, seed, keyargs)
        push!(seeds, seed)
        visited[seed] = length(seeds)
        # new seed thresh
        seed_weights = weights[getedgeindex.(seed, 1:3)]
        new_seed_thresh = NBINS - div(NBINS - sum(seed_weights)/3, 2)
        return new_seed_thresh
    end
    return addseed!
end

function getseedqueue(weights)
    queue = PQueue{Int}(3NBINS)
    for (i, w) in enumerate(sum([w == 0 ? UInt8(255) : w for w in weights]; dims=1))
        if w != 255
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

function growRegionUnwrap!(wrapped, weights, nbins, keyargs, max_seeds=2)
    stridelist = strides(wrapped)
    visited = zeros(UInt8, size(wrapped))
    notvisited(i) = checkbounds(Bool, visited, i) && (visited[i] == 0)
    pqueue = PQueue{Int}(nbins)
    seeds = Int[]
    seedqueue = getseedqueue(weights, nbins)

    function addneighbors(newvox)
        for i in 1:6 # 6 directions
            e = getnewedge(newvox, notvisited, stridelist, i)
            if e != 0 && weights[e] > 0
                push!(pqueue, e, weights[e])
            end
        end
    end


    function addseed!(seeds, seedqueue, pqueue, visited, weights, nbins)
        seed = findseed!(seedqueue, weights, visited)
        if seed == 0
            return 255
        end
        seedcorrection!(wrapped, seed, keyargs)
        push!(seeds, seed)
        addneighbors(seed)
        visited[seed] = length(seeds)
        weight_thresh = nbins - div(nbins - weights[seed], 2)
        return weight_thresh
    end
    weight_thresh = addseed!(seeds, seedqueue, pqueue, visited, weights, nbins)


    #@show Int.(weights)
    while !isempty(pqueue)
        if length(seeds) < max_seeds && pqueue.min > weight_thresh
            weight_thresh = addseed!(seeds, seedqueue, pqueue, visited, weights, nbins)
        end
        edge = pop!(pqueue)
        oldvox, newvox = getvoxelsfromedge(edge, visited, stridelist)
        #@show edge oldvox newvox
        if visited[newvox] == 0
            unwrapedge!(wrapped, oldvox, newvox)
            visited[newvox] = visited[oldvox]
            addneighbors(newvox)
        end
    end

    correct_regions!(wrapped, visited, length(seeds))
    return wrapped
end

function correct_regions!(wrapped, visited, nregions)
    offsets = zeros(nregions, nregions)
    offset_counts = zeros(Int, nregions, nregions)
    stridelist = strides(wrapped)
    for dim in 1:3
        neighbor = stridelist[dim]
        for I in LinearIndices(wrapped)
            J = I + neighbor
            offsets[visited[I], visited[J]] += wrapped[I] - wrapped[J]
            offset_counts[visited[I], visited[J]] += 1
        end
    end
    for i in 1:nregions, j in i:nregions
        offset_counts[i,j] += offset_counts[j,i]
        offset_counts[j,i] = 0
        offsets[i,j] -= offsets[j,i]
        offsets[j,i] = -offsets[i,j]
    end
    corrected = falses(nregions)
    corrected[1] = true
    while !all(corrected)
        (_, (i,j)) = findmax(offset_counts) # most connections # i<j
        if !corrected[j] && offset_counts != 0
            offset = round((offsets[i,j] / offset_counts) / 2π)
            if offset != 0
                wrapped[visited .== j] .-= offset * 2π
                visited[visited .== j] .= i
            end
            corrected[j] = true
        end
        offset_counts[i,j] = 0
    end
end

function addseed!(seeds, seedqueue, pqueue, visited, weights, nbins, getnewseed!)
    seed = getnewseed!(seedqueue, visited)
    if seed == 0
        return 255
    end
    push!(seeds, seed)
    addneighbors(seed)
    visited[seed] = length(seeds)
    weight_thresh = nbins - div(nbins - weights[seed], 2)
    return weight_thresh
end

function getvoxelsfromedge(edge, visited, stridelist)
    dim = getdimfromedge(edge)
    vox = getfirstvoxfromedge(edge)
    neighbor = vox + stridelist[dim] # direct neigbor in dim
    if visited[neighbor] == 0
        return vox, neighbor
    else
        return neighbor, vox
    end
end

# edge calculations
getdimfromedge(edge) = (edge - 1) % 3 + 1
getfirstvoxfromedge(edge) = div(edge - 1, 3) + 1
getedgeindex(leftvoxel, dim) = dim + 3(leftvoxel-1)

function unwrapedge!(wrapped, oldvox, newvox)
    wrapped[newvox] = unwrapvoxel(wrapped[newvox], wrapped[oldvox])
end
unwrapvoxel(new, old) = new - 2pi * round((new - old) / 2pi)

function getnewedge(v, notvisited, stridelist, i)
    iDim = div(i+1,2)
    n = stridelist[iDim] # neigbor-offset in dimension iDim
    if iseven(i)
        if notvisited(v+n) getedgeindex(v, iDim) else 0 end
    else
        if notvisited(v-n) getedgeindex(v-n, iDim) else 0 end
    end
end

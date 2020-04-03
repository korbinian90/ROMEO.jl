function growRegionUnwrap!(wrapped, weights, nbins, keyargs, maxseeds=50)
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
        seed_weights = weights[getedgeindex.(seed, 1:3)]
        weight_thresh = nbins - div(nbins - maximum(seed_weights), 2)
        #@show Int.(seed_weights) weight_thresh sum(visited .== 0)
        return weight_thresh
    end
    weight_thresh = addseed!(seeds, seedqueue, pqueue, visited, weights, nbins)

    #@show Int.(weights)
    while !isempty(pqueue)
        if length(seeds) < maxseeds && pqueue.min > weight_thresh
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
    #@show length(seeds)
    if !(haskey(keyargs, :phase2) && haskey(keyargs, :TEs))
        @show "correct regions"
        correct_regions!(wrapped, copy(visited), length(seeds), weights)
    end
    #return wrapped
    return wrapped, visited, weights
end

function correct_regions!(wrapped, visited, nregions, weights)
    # biggest_region = findmax(countmap(visited; algo=:dict))[2]
    offsets = zeros(nregions, nregions)
    offset_counts = zeros(Int, nregions, nregions)
    stridelist = strides(wrapped)
    for dim in 1:3
        neighbor = CartesianIndex(ntuple(i->i==dim ? 1 : 0, 3))
        for I in CartesianIndices(wrapped)
            J = I + neighbor
            if checkbounds(Bool, wrapped, J)# && weights[dim, I] < 255
                ri = visited[I]
                rj = visited[J]
                if ri != 0 && rj != 0
                    offsets[ri, rj] += wrapped[I] - wrapped[J]
                    offset_counts[ri, rj] += 1
                end
            end
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
    while any(offset_counts .> 0)
        (_, I) = findmax(offset_counts) # most connections
        (i,j) = Tuple(I) # i<j
        if !corrected[j] && offset_counts[i,j] > 0
            offset = round((offsets[i,j] / offset_counts[i,j]) / 2π)
            if offset != 0
                println("region $j adjusted to $i with offset $offset * 2π")
                println("offset_count $(offset_counts[i,j])")
                wrapped[visited .== j] .+= offset * 2π
                visited[visited .== j] .= i
            end
            corrected[j] = true
        end
        offset_counts[i,j] = -1
    end
    #=
    for region in unique(visited)
        offset =
    end
    =#
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

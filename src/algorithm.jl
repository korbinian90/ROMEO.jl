function growRegionUnwrap!(wrapped, weights, seed, nbins)
    stridelist = strides(wrapped)
    visited = falses(size(wrapped))
    pqueue = PQueue(nbins, seed)

    while !isempty(pqueue)
        edge = pop!(pqueue)
        oldvox, newvox = getvoxelsfromedge(edge, visited, stridelist)
        if !visited[newvox]
            unwrapedge!(wrapped, oldvox, newvox)
            visited[newvox] = true
            for e in getnewedges(newvox, visited, stridelist)
                if weights[e] > 0
                    push!(pqueue, e, weights[e])
                end
            end
        end
    end
    return wrapped
end

function getvoxelsfromedge(edge, visited, stridelist)
    dim = getdimfromedge(edge)
    vox = getfirstvoxfromedge(edge)
    neighbor = vox + stridelist[dim] # direct neigbor in dim
    if !visited[neighbor]
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

function getnewedges(v, visited, stridelist)
    notvisited(i) = checkbounds(Bool, visited, i) && !visited[i]

    edges = []
    for iDim = 1:3
        n = stridelist[iDim] # neigbor-offset in dimension iDim
        if notvisited(v+n) push!(edges, getedgeindex(v, iDim)) end
        if notvisited(v-n) push!(edges, getedgeindex(v-n, iDim)) end
    end
    return edges
end

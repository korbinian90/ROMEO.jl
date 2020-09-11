function grow_region_unwrap!(
    wrapped, weights, visited=zeros(UInt8, size(wrapped)), pqueue=PQueue{Int}(NBINS);
    maxseeds=1, merge_regions=false, correct_regions=false, wrap_addition=0, keyargs...
    )
    ## Init
    maxseeds = min(255, maxseeds) # Stored in UInt8
    dimoffsets = getdimoffsets(wrapped)
    notvisited(i) = checkbounds(Bool, visited, i) && (visited[i] == 0)
    seeds = Int[]
    new_seed_thresh = 256
    if isempty(pqueue) # no seed added yet
        addseed! = getseedfunction(seeds, pqueue, visited, weights, wrapped, keyargs)
        new_seed_thresh = addseed!()
    end
    ## MST loop
    while !isempty(pqueue)
        if length(seeds) < maxseeds && pqueue.min > new_seed_thresh
            new_seed_thresh = addseed!()
        end
        edge = dequeue!(pqueue)
        oldvox, newvox = getvoxelsfromedge(edge, visited, dimoffsets)
        if visited[newvox] == 0
            unwrapedge!(wrapped, oldvox, newvox, visited, wrap_addition)
            visited[newvox] = visited[oldvox]
            for i in 1:6 # 6 directions
                e = getnewedge(newvox, notvisited, dimoffsets, i)
                if e != 0 && weights[e] > 0
                    enqueue!(pqueue, e, weights[e])
                end
            end
        end
    end
    ## Region merging
    if merge_regions regions = merge_regions!(wrapped, visited, length(seeds), weights) end
    if merge_regions && correct_regions correct_regions!(wrapped, visited, regions) end
    return visited
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

function unwrapedge!(wrapped, oldvox, newvox, visited, x)
    oo = 2oldvox - newvox
    d = 0
    if checkbounds(Bool, wrapped, oo) && visited[oo] != 0 # neighbor behind is visited
        v = wrapped[oldvox] - wrapped[oo]
        d = if v < -x # threshold
            -x
        elseif v > x
            x
        else
            v
        end
    end
    wrapped[newvox] = unwrapvoxel(wrapped[newvox], wrapped[oldvox] + d)
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

function grow_region_unwrap!(
    wrapped, weights, visited=zeros(UInt8, size(wrapped)), pqueue=PQueue{Int}(NBINS);
    maxseeds=1, merge_regions=false, correct_regions=false, wrap_addition=0, keyargs...
    )
    ## Init
    stridelist = strides(wrapped)
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
        edge = pop!(pqueue)
        oldvox, newvox = getvoxelsfromedge(edge, visited, stridelist)
        if visited[newvox] == 0
            unwrapedge!(wrapped, oldvox, newvox, visited, wrap_addition)
            visited[newvox] = visited[oldvox]
            for i in 1:6 # 6 directions
                e = getnewedge(newvox, notvisited, stridelist, i)
                if e != 0 && weights[e] > 0
                    push!(pqueue, e, weights[e])
                end
            end
        end
    end
    ## Region merging
    if merge_regions regions = merge_regions!(wrapped, visited, length(seeds), weights) end
    if correct_regions correct_regions!(wrapped, visited, regions) end
end

function getseedfunction(seeds, pqueue, visited, weights, wrapped, keyargs)
    seedqueue = getseedqueue(weights)
    notvisited(i) = checkbounds(Bool, visited, i) && (visited[i] == 0)
    stridelist = strides(wrapped)
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

function correct_regions!(wrapped, visited, regions)
    for r in regions
        wrapped[visited .== r] .-= (2π * median(round.(wrapped[visited .== r] ./ 2π)))
    end
end

function merge_regions!(wrapped, visited, nregions, weights)
    mask = sum(weights; dims=1)
    region_size = countmap(visited) # TODO could use weight instead
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
                if ri != 0 && rj != 0 && ri != rj
                    w = 255 - weights[dim,I]
                    if w == 255 w = 0 end
                    offsets[ri, rj] += (wrapped[I] - wrapped[J]) * w
                    offset_counts[ri, rj] += w
                end
            end
        end
    end
    for i in 1:nregions, j in i:nregions
        offset_counts[i,j] += offset_counts[j,i]
        offset_counts[j,i] = offset_counts[i,j]
        offsets[i,j] -= offsets[j,i]
        offsets[j,i] = -offsets[i,j]
    end
    corrected = falses(nregions)
    remaining_regions = Int[]
    while !all(corrected)
        largest_uncorrected_region = try
            findmax(filter(p -> first(p) != 0 && !corrected[first(p)], region_size))[2]
        catch
            @show region_size size(corrected) nregions
            throw(error())
        end
        # TODO correct region? No, offsets are already calculated
        corrected[largest_uncorrected_region] = true
        push!(remaining_regions, largest_uncorrected_region)

        # TODO multiple rounds until no change?
        sorted_offsets = get_offset_count_sorted(offset_counts, corrected)
        for I in sorted_offsets
            (i,j) = Tuple(I)
            offset = round((offsets[I] / offset_counts[I]) / 2π)
            if offset != 0
                wrapped[visited .== j] .+= offset * 2π
                visited[visited .== j] .= i
                # TODO region merging offset_count and offset calculation
            end
            corrected[j] = true
            offset_counts[i,j] = offset_counts[j,i] = -1
        end
    end
    return remaining_regions
end

function get_offset_count_sorted(offset_counts, corrected)
    f(I) = corrected[I[1]] && !corrected[I[2]]
    sort(filter(f, findall(offset_counts .> 0)), by=i->offset_counts[i], rev=true)
end
initqueue(seed::Int, weights) = initqueue([seed], weights)
function initqueue(seeds, weights)
    pq = PQueue{eltype(seeds)}(NBINS)
    for seed in seeds
        push!(pq, seed, weights[seed])
    end
    return pq
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

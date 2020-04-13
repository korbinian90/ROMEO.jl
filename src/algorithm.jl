function growRegionUnwrap!(wrapped, weights, seed, nbins, data)
    stridelist = strides(wrapped)
    visited = zeros(UInt8, size(wrapped))
    notvisited(i) = checkbounds(Bool, visited, i) && (visited[i] == 0)
    pqueue = PQueue(nbins, seed)
    visited[getfirstvoxfromedge(seed)] = 1
    #@show typeof(wrapped) typeof(data) typeof(stridelist)
    while !isempty(pqueue)
        edge = pop!(pqueue)
        oldvox::Int, newvox::Int = getvoxelsfromedge(edge, visited, stridelist)
        if visited[newvox] == 0
            unwrapedge!(wrapped, oldvox::Int, newvox::Int, data, stridelist)
            visited[newvox] = visited[oldvox]
            for i in 1:6 # 6 directions
                e = getnewedge(newvox, notvisited, stridelist, i)
                if e != 0 && weights[e] > 0
                    push!(pqueue, e, weights[e])
                end
            end
        end
    end
    return wrapped
end

function getnewseed!(wrapped, weights, visited; keyargs...)
    seed = findseed(wrapped, weights, visited)
    if haskey(keyargs, :phase2) && haskey(keyargs, :TEs) # requires multiecho
        seedcorrection!(wrapped, seed, keyargs[:phase2], keyargs[:TEs])
    end
    visited
    return seed
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

function unwrapedge!(wrapped, oldvox::Int, newvox::Int, d::Nothing, stridelist)
    wrapped[newvox] = unwrapvoxel(wrapped[newvox], wrapped[oldvox])
end

function unwrapedge!(wrapped, oldvox::Int, newvox::Int, d::Data, stridelist::NTuple{4,Int64})::Nothing
    wrapped[newvox] = unwrapvoxel(wrapped[newvox], wrapped[oldvox])
    TE = d.TEs
    neco = size(wrapped, 4)
    n = stridelist[4] # neigbor-offset in dimension iDim
    echoinds = (0:neco-1) .* n
    oldvox4D = oldvox .+ echoinds
    newvox4D = newvox .+ echoinds
    if ndims(d.mag) == 3
        echoinds = zeros(Int, neco)
    end
    oldvox4D_mag = oldvox .+ echoinds
    newvox4D_mag = newvox .+ echoinds

    #PO = d.PO[oldvox]
    #B0 = d.B0[oldvox]
    #newB0 = wrapped[newvox] / TE[1]
    #newB0weight = d.mag[newvox] * TE[1]
    est = (wrapped[newvox] - wrapped[oldvox]) / TE[1]
    est_weight = d.mag[newvox] * TE[1]

    for i in 2:neco
        #val = median((wrapped[oldvox4D[i]], PO + TE[i] * B0, PO + TE[i] * newB0))
        val = wrapped[oldvox4D[i]] + est * TE[i]
        wrapped[newvox4D[i]] = unwrapvoxel(wrapped[newvox4D[i]], val)

        diff = wrapped[newvox4D[i]] - wrapped[oldvox4D[i]]
        diff_weighted = diff * d.mag[newvox4D_mag[i]]
        w_sum = est_weight + d.mag[newvox4D_mag[i]] * TE[i]

        est = (est * est_weight + diff_weighted) / w_sum
        est_weight += d.mag[newvox4D_mag[i]] * TE[i]
        #wi = d.mag[newvox4D[i]] * TE[i]
        #newB0 = (newB0 * newB0weight + wrapped[newvox4D[i]] / TE[i] * wi) / (newB0weight + wi)
        #newB0weight = newB0weight + wi
    end

    #d.B0[newvox] = newB0
    #d.PO[newvox] = median(wrapped[newvox4D] .- newB0 .* TE)
    return nothing
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

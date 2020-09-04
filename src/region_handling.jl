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
    stridelist = getdimoffsets(wrapped)
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
    return sort(filter(f, findall(offset_counts .> 0)), by=i->offset_counts[i], rev=true)
end

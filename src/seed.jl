function findseed(wrapped, weights)
    cp = copy(weights)
    cp[cp .== 0] .= 255
    filtered = dilate(cp, 2:4)
    (_, ind) = findmin(filtered)
    return LinearIndices(weights)[ind]
end

function findseed(wrapped, weights, visited)
    cp = copy(weights)
    cp[(visited .!== 0) .& (cp .== 0)] .= 255
    (_, ind) = findmin(cp)
    return LinearIndices(weights)[ind]
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

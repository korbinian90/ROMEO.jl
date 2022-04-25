function wrapfit(phase, TEs, num=4)
    best_fiterr = Inf
    best_unwrapped = phase
    for off in -num*2pi:2pi:num*2pi
        @show off
        @show unwrapped = fit(phase, off, TEs)
        @show fiterr = fiterror(unwrapped, TEs) + sqrt(abs(off)) / 100
        if fiterr < best_fiterr
            best_fiterr = fiterr
            best_unwrapped = unwrapped
        end
    end
    return best_unwrapped
end

function fit(phase, off, TEs)
    unwrapped = copy(phase)
    b0 = (unwrapped[1] + off) / TEs[1]
    for i in 1:length(TEs)
        unwrapped[i] = unwrapvoxel(unwrapped[i], b0 * TEs[i])
        b0 = (b0 + unwrapped[i] / TEs[i]) / 2
    end
    return unwrapped
end

function fiterror(phase, TEs)
    b0 = phase ./ TEs
    mean_b0 = sum(b0) / length(b0)
    error = sum((b0 .- mean_b0) .^2)
    return error
end

function detect_region_for_echo_correction(phase; keyargs...)
    kw = Dict()
    if haskey(keyargs, :mag)
        kw[:mag] = keyargs[:mag]
        if size(kw[:mag],4) > 1
            kw[:mag] = kw[:mag][:,:,:,end]
        end
    end
    kw[:phase2] = phase[:,:,:,end-1]
    kw[:TEs] = keyargs[:TEs][end:-1:end-1]
    if haskey(keyargs, :mask)
        kw[:mask] = keyargs[:mask]
    end

    qm = calculateweights(phase[:,:,:,end]; kw...)
    qm = mapwindow(minimum, qm, (5,5,3))
    mask = qm .> quantile(qm, 0.8)
    med = median(phase[:,:,:,1][mask])
    locations = phase .== med .& mask
    @show vox = findfirst(locations)
    return vox
end

function echo_jump_correction(phase; keyargs...)
    vox = detect_region_for_echo_correction(phase; keyargs...)
    # TODO area around vox (3x3x3 cube?)
    uw = wrapfit(phase[vox...,:], keyargs[:TEs])
    offset = ones(1,1,1,size(phase,4))
    offset .= uw .- phase[vox...,:]
    return phase .+ offset
end

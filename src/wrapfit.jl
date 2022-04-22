function wrapfit(phase, TEs, num=4)
    best_fiterr = Inf
    best_unwrapped = phase
    for off in -num*2pi:2pi:num*2pi
        @show off
        @show unwrapped = fit(phase, off, TEs)
        @show fiterr = fiterror(unwrapped, TEs) + abs(off) / 100
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

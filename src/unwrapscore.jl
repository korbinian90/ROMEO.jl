function unwrapscore(phase, mag=ones(size(phase)), λ=0.01)
    score = 0
    dims = ndims(phase)
    sz = size(phase)
    for idim in 1:dims
        s0 = getshiftedindices(sz, idim, 0)
        s1 = getshiftedindices(sz, idim, 1)
        s2 = getshiftedindices(sz, idim, 2)
        score += sum(
            abs.(phase[s0...] .- 2phase[s1...] .+ phase[s2...]) .*
            (mag[s0...] .+ mag[s1...] .+ mag[s2...])
        ) / length(phase)
    end
    score += λ * sqrt(sum(phase.^2) / length(phase))
    return score
end

# size reduced by two
function getshiftedindices(sz, dim, shift)
    return (fill(:, dim-1)..., shift .+ (1:sz[dim]-2), fill(:, length(sz)-dim)...)
end

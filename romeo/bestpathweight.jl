function getD(φ)
    directions = Iterators.product(-1:1,-1:1,-1:1)
    neigbors = unique(abs(sum(strides(φ) .* d)) for d in directions if d != (0,0,0))
    D2 = zeros(size(φ))
    @inbounds for n in neigbors, i in 1+n:length(φ)-n
        D2[i] += (γ(φ[i-n] - φ[i]) - γ(φ[i] - φ[i+n])) ^ 2
    end
    return sqrt.(D2)
end

function γ(x) # fastest version
    if x < -π
        x+typeof(x)(2π)
    elseif x > π
        x-typeof(x)(2π)
    else
        x
    end
end

function getbestpathweight(φ)
    R = 1 ./ getD(φ)
    weight = zeros(3, size(R)...)
    for idim in 1:3
        n = strides(R)[idim]
        for i in 1:length(R)-n
            weight[idim + 3i] = R[i] + R[i+n]
        end
    end
    return weight
end

function γ(x::AbstractFloat) # faster if only one wrap can occur
    if x < -π
        x+typeof(x)(2π)
    elseif x > π
        x-typeof(x)(2π)
    else
        x
    end
end

getdimoffsets(A) = (1, cumprod(size(A)[1:end-1])...)

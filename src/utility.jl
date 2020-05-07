function γ(x::AbstractFloat) # faster if only one wrap can occur
    if x < -π
        x+typeof(x)(2π)
    elseif x > π
        x-typeof(x)(2π)
    else
        x
    end
end

# makes strides available for BitArray and NIVolume
Base.strides(A::AbstractArray) = (1, cumprod(collect(size(A)))[1:end-1]...)

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

# Statistics.median! without the repr of the input in its error message, which
# static compilation cannot resolve. Same algorithm, so results are identical.
function _median!(v::AbstractVector)
    isempty(v) && throw(ArgumentError("median of an empty array is undefined"))
    nanix = findfirst(isnan, v)
    isnothing(nanix) || return v[nanix]
    n = length(v)
    mid = div(1 + n, 2)
    if isodd(n)
        return Statistics.middle(partialsort!(v, mid))
    else
        m = partialsort!(v, mid:mid+1)
        return Statistics.middle(m[1], m[2])
    end
end

# Keyword arguments for one echo of a multi-echo array: the reference echo as
# phase2 with both echo times, and the magnitude of that echo. A NamedTuple
# instead of a Dict{Symbol,Any} keeps every call statically typed.
function echo_keyargs(keyargs, ieco, phase2, TEs)
    args = merge(values(keyargs), (; phase2, TEs))
    if haskey(args, :mag)
        args = merge(args, (; mag=select_echo(args.mag, ieco)))
    end
    return args
end
select_echo(mag::AbstractArray{<:Any,4}, ieco) = mag[:,:,:,ieco]
select_echo(mag, ieco) = mag

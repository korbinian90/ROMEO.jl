function unwrap!(wrapped::AbstractArray{T,3}; regions=zeros(UInt8, size(wrapped)), keyargs...) where T
    weights = calculateweights(wrapped; keyargs...)
    @assert sum(weights) != 0 "Unwrap-weights are all zero!"

    regions .= grow_region_unwrap!(wrapped, weights; keyargs...)

    if haskey(keyargs, :correctglobal) && keyargs[:correctglobal]
        mask = if haskey(keyargs, :mask)
            keyargs[:mask]
        else
            trues(size(wrapped))
        end
        wrapped .-= (2π * median(round.(filter(isfinite, wrapped[mask]) ./ 2π))) # TODO time (sample)
    end
    return wrapped
end

function unwrap!(wrapped::AbstractArray; keyargs...)
    sz = size(wrapped)
    @assert ndims(wrapped) <= 3 "This is 3D (4D) unwrapping! data is $(ndims(wrapped))D"
    if ndims(wrapped) <= 2 # algorithm requires 3D input
        wrapped = reshape(wrapped, size(wrapped)..., ones(Int, 3-ndims(wrapped))...)
    end
    unwrap!(wrapped; keyargs...)
    return reshape(wrapped, sz)
end

"""
    unwrap(wrapped::AbstractArray; keyargs...)

ROMEO unwrapping for 3D and 4D data.

###  Optional keyword arguments:

- `TEs`: Required for 4D data. The echo times for multi-echo data. In the case of single-echo 
    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).
- `weights`: Options are [`:romeo`] | `:romeo2` | `:romeo3` | `:bestpath`.
- `mag`: The magnitude is used to improve the unwrapping-path.
- `mask`: Unwrapping is only performed inside the mask.
- `phase2`: A second reference phase image (possibly with different echo time).
    It is used for calculating the phasecoherence weight. This is automatically
    done for 4D multi-echo input and therefore not required.
- `correctglobal=false`: If `true` corrects global n2π offsets.
- `individual=false`: If `true` perform individual unwrapping of echos.
    Type `?unwrap_individual` for more information
- `template=1`: echo that is spatially unwrapped (if `individual` is `false`)
- `maxseeds=1`: higher values allow more seperate regions
- `merge_regions=false`: spatially merge neighboring regions after unwrapping
- `correct_regions=false`: bring each regions median closest to 0 by adding n2π
- `wrap_addition=0`: [0;π], allows 'linear unwrapping', neighbors can have more
    (π+wrap_addition) phase difference
- `temporal_uncertain_unwrapping=false`: uses spatial unwrapping on voxels that
    have high uncertainty values after temporal unwrapping

# Examples
```julia-repl
julia> using MriResearchTools
julia> phase = readphase("phase_3echo.nii")
julia> unwrapped = unwrap(phase; TEs=[1,2,3])
julia> savenii(unwrapped, "unwrapped.nii"; header=header(phase))
```
"""
unwrap, unwrap!

unwrap(wrapped; keyargs...) = unwrap!(copy(wrapped); keyargs...)

function unwrap!(wrapped::AbstractArray{T,4}; TEs, individual=false,
        template=1, p2ref=ifelse(template==1, 2, template-1),
        temporal_uncertain_unwrapping=false, keyargs...) where T
    if individual return unwrap_individual!(wrapped; TEs, keyargs...) end
    ## INIT
    args = Dict{Symbol, Any}(keyargs)
    args[:phase2] = wrapped[:,:,:,p2ref]
    args[:TEs] = TEs[[template, p2ref]]
    if haskey(args, :mag)
        args[:mag] = args[:mag][:,:,:,template]
    end
    ## Calculate
    weights = calculateweights(view(wrapped,:,:,:,template); args...)
    unwrap!(view(wrapped,:,:,:,template); args..., weights) # rightmost keyarg takes precedence
    quality = similar(wrapped)
    V = falses(size(wrapped))
    for ieco in [(template-1):-1:1; (template+1):length(TEs)]
        iref = if (ieco < template) ieco+1 else ieco-1 end
        refvalue = wrapped[:,:,:,iref] .* (TEs[ieco] / TEs[iref])
        w = view(wrapped,:,:,:,ieco)
        w .= unwrapvoxel.(w, refvalue) # temporal unwrapping
        
        if temporal_uncertain_unwrapping # TODO extract as function
            quality[:,:,:,ieco] .= getquality.(w, refvalue)
            visited = quality[:,:,:,ieco] .< π/2
            mask = if haskey(keyargs, :mask)
                keyargs[:mask]
            else
                dropdims(sum(weights; dims=1); dims=1) .< 100
            end
            visited[.!mask] .= true
            V[:,:,:,ieco] = visited
            if any(visited) && !all(visited)
                edges = getseededges(visited)
                edges = filter(e -> weights[e] != 0, edges)
                grow_region_unwrap!(w, weights, visited, initqueue(edges, weights))
            end
        end
    end
    return wrapped#, quality, weights, V
end

function getquality(vox, ref)
    return abs(vox - ref)
end

function getseededges(visited::BitArray)
    stridelist = getdimoffsets(visited)
    edges = Int64[]
    for dim in 1:3, I in LinearIndices(visited)
        J = I + stridelist[dim]
        if checkbounds(Bool, visited, J) # borders should be no problem due to invalid weights
            if visited[I] + visited[J] == 1 # one visited and one not visited
                push!(edges, getedgeindex(I, dim))
            end
        end
    end
    return edges
end

initqueue(seed::Int, weights) = initqueue([seed], weights)
function initqueue(seeds, weights)
    pq = PQueue{eltype(seeds)}(NBINS)
    for seed in seeds
        enqueue!(pq, seed, weights[seed])
    end
    return pq
end

"""
    unwrap_individual(wrapped::AbstractArray{T,4}; TEs, keyargs...) where T

Performs individual unwrapping of the echoes instead of temporal unwrapping.
Still uses multi-echo information to improve the quality map.
This function is identical to `unwrap` with the flag `individual=true`.
The syntax is identical to unwrap, but doesn't support the `temporal_uncertain_unwrapping` and `template` options:
$(@doc unwrap)

"""
unwrap_individual, unwrap_individual!

unwrap_individual(wrapped; keyargs...) = unwrap_individual!(copy(wrapped); keyargs...)
function unwrap_individual!(wrapped::AbstractArray{T,4}; TEs, keyargs...) where T
    args = Dict{Symbol,Any}(keyargs)
    Threads.@threads for i in 1:length(TEs)
        e2 = if (i == 1) 2 else i-1 end
        if haskey(keyargs, :mag) args[:mag] = keyargs[:mag][:,:,:,i] end
        unwrap!(view(wrapped,:,:,:,i); phase2=wrapped[:,:,:,e2], TEs=TEs[[i,e2]], args...)
    end
    if haskey(keyargs, :correctglobal) && keyargs[:correctglobal]
        correct_multi_echo_wraps!(wrapped; TEs, keyargs...)
    end
    return wrapped
end

function correct_multi_echo_wraps!(wrapped; TEs, mask=trues(size(wrapped)), keyargs...)
    for ieco in 2:length(TEs)
        iref = ieco - 1
        nwraps = median(round.((filter(isfinite, wrapped[:,:,:,iref][mask]) .* (TEs[ieco] / TEs[iref]) .- filter(isfinite, wrapped[:,:,:,ieco][mask])) / 2π))
        wrapped[:,:,:,ieco] .+= 2π * nwraps 
    end
end

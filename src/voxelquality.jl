"""
    voxelquality(phase::AbstractArray; keyargs...)

Calculates a quality for each voxel. The voxel quality can be used to create a mask.

# Examples
```julia-repl
julia> qmap = voxelquality(phase_3echo; TEs=[1,2,3]);
julia> mask = robustmask(qmap);
```
     
Takes the same inputs as `romeo`:
$(@doc romeo)

See also [`romeo`](@ref)
""" 
function voxelquality(phase; keyargs...)
    weights = calculateweights(phase; type=Float32, rescale=x->x, keyargs...) # [0;1]
    qmap = dropdims(sum(weights; dims=1); dims=1)
    qmap[2:end,:,:] .+= weights[1,1:end-1,:,:]
    qmap[:,2:end,:] .+= weights[2,:,1:end-1,:]
    qmap[:,:,2:end] .+= weights[3,:,:,1:end-1]
    return qmap ./ 6 # [0;1]
end

function calculateweights(phase::AbstractArray{T,4}; TEs, template=1, p2ref=2, keyargs...) where T
    args = Dict{Symbol, Any}(keyargs)
    args[:phase2] = phase[:,:,:,p2ref]
    args[:TEs] = TEs[[template, p2ref]]
    if haskey(args, :mag) && size(args[:mag], 4) > 1
        args[:mag] = args[:mag][:,:,:,template]
    end
    return calculateweights(view(phase,:,:,:,template); args...)
end

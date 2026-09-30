"""
    voxelquality(phase::AbstractArray; keyargs...)

Calculates a quality for each voxel. The voxel quality can be used to create a mask.
The quality range is [0;1]

# Examples
```julia-repl
julia> qmap = voxelquality(phase_3echo; TEs=[1,2,3]);
julia> mask = robustmask(qmap);
```
     
Takes the same inputs as `romeo`/`unwrap`:

See also [`unwrap`](@ref)
""" 
function voxelquality(phase; keyargs...)
    # convert: a bestpath or user-supplied weights array keeps its own element type
    weights = convert(Array{Float32,4}, calculateweights(phase; type=Val(Float32), rescale=identity, keyargs...)) # [0;1]
    qmap = dropdims(sum(weights; dims=1); dims=1)
    qmap[2:end,:,:] .+= weights[1,1:end-1,:,:]
    qmap[:,2:end,:] .+= weights[2,:,1:end-1,:]
    qmap[:,:,2:end] .+= weights[3,:,:,1:end-1]
    return qmap ./ 6 # [0;1]
end

function calculateweights(phase::AbstractArray{T,4}; TEs, template=1, p2ref=2, keyargs...) where T
    args = echo_keyargs(keyargs, template, phase[:,:,:,p2ref], TEs[[template, p2ref]])
    return calculateweights(view(phase,:,:,:,template); args...)
end

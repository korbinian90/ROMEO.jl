## Phase singularity (residue / vortex line) detection
#
# A residue is a non-zero circulation of the wrapped phase differences around a
# 2x2 plaquette. Where it occurs, no integer field n(x) exists that makes
# φ + 2πn continuous everywhere. This is a topological obstruction and not a
# shortcoming of the unwrapping algorithm. Physically a singularity sits where
# the complex signal vanishes (Re(S) = 0 and Im(S) = 0), which is a line in 3D.
# The dominant cause in high resolution data is partial voluming at a vessel
# wall, where two compartments with ~π phase difference cancel each other.
#
# In 3D the residue field is divergence free (the charges of the 6 faces of a
# voxel cube sum to zero), therefore residues form closed lines, around a vein
# typically rings that enclose the vessel.
#
# The residues only depend on the wrapped phase differences and are therefore
# invariant under adding n2π to individual voxels. They can be computed from
# the wrapped input as well as from the congruent ROMEO output.

# normal dimension -> the two dimensions spanned by the plaquette (cyclic, so
# that the residue field is divergence free)
const PLAQUETTE_DIMS = ((2, 3), (3, 1), (1, 2))

unitindex(d) = CartesianIndex(ntuple(i -> Int(i == d), 3))

tomask(mask) = isnothing(mask) ? nothing : BitArray(mask .!= 0)
tomask(mask::BitArray) = mask

# 1D and 2D input is handled as 3D with singleton dimensions, like `unwrap`
const LowDim{T} = Union{AbstractArray{T,1},AbstractArray{T,2}}
to3d(x) = reshape(x, size(x)..., ntuple(_ -> 1, 3 - ndims(x))...)
to3d(::Nothing) = nothing

"""
    residues(phase::AbstractArray{<:Real,3}; mask=nothing)

Charge of every 2x2 plaquette of the phase (Goldstein residues, extended to 3D).

`residues(phase)[c,i,j,k]` is the charge (`-1`, `0` or `1`) of the plaquette
with the corner voxel `(i,j,k)` that is normal to the dimension `c`. The
returned field is divergence free, non-zero plaquettes therefore form closed
lines (see [`detect_singularities`](@ref)).

Wrapped and unwrapped (congruent) phase give the identical result.
"""
function residues(phase::AbstractArray{<:Real,3}; mask=nothing)
    sz = size(phase)
    mask = tomask(mask)
    res = zeros(Int8, 3, sz...)
    Threads.@threads for k in 1:sz[3] # threading over slices to avoid write conflicts
        for c in 1:3
            a, b = PLAQUETTE_DIMS[c]
            (k == sz[3] && (a == 3 || b == 3)) && continue
            ea, eb = unitindex(a), unitindex(b)
            ranges = (1:(sz[1] - (a == 1 || b == 1)), 1:(sz[2] - (a == 2 || b == 2)), k:k)
            for I in CartesianIndices(ranges)
                if !isnothing(mask) && !(mask[I] && mask[I+ea] && mask[I+eb] && mask[I+ea+eb])
                    continue
                end
                s = rem2pi(phase[I+ea] - phase[I], RoundNearest) +
                    rem2pi(phase[I+ea+eb] - phase[I+ea], RoundNearest) +
                    rem2pi(phase[I+eb] - phase[I+ea+eb], RoundNearest) +
                    rem2pi(phase[I] - phase[I+eb], RoundNearest)
                if isfinite(s)
                    res[c, I] = round(Int8, s / 2π)
                end
            end
        end
    end
    return res
end

"""
    SingularityLoop

One connected component of the residue field, i.e. one vortex line.

- `plaquettes`: linear indices into the `(3, size(phase)...)` residue array
- `voxels`: linear indices of the voxels touched by these plaquettes
- `closed`: `true` if the line closes inside the volume (and mask). Open lines
    end at the volume or mask border and cannot be corrected consistently.
"""
struct SingularityLoop
    plaquettes::Vector{Int}
    voxels::Vector{Int}
    closed::Bool
end

loop_length(l::SingularityLoop) = length(l.plaquettes)

"""
    Singularities

Result of [`detect_singularities`](@ref): the residue field and the vortex
lines (`loops`) it forms.
"""
struct Singularities
    residues::Array{Int8,4}
    loops::Vector{SingularityLoop}
    size::NTuple{3,Int}
end

"""
    detect_singularities(phase::AbstractArray{<:Real,3}; mask=nothing)

Detect phase singularities and group them into vortex lines.

The phase can be the wrapped input or the unwrapped ROMEO output, both give the
same result. Detection is parameter free and O(N), it is the definition of the
problem rather than a heuristic.

# Examples
```julia-repl
julia> unwrapped = unwrap(phase; mag);
julia> s = detect_singularities(unwrapped);
julia> sum(abs, s.residues) # number of residue plaquettes
julia> savenii(singularity_mask(s), "singularities.nii"; header=header(phase));
```

See also [`singularity_mask`](@ref), [`branchcuts`](@ref),
[`fix_singularities!`](@ref)
"""
function detect_singularities(phase::AbstractArray{<:Real,3}; mask=nothing)
    res = residues(phase; mask)
    return Singularities(res, find_loops(res), size(phase))
end

detect_singularities(phase::AbstractArray{<:Real,4}; keyargs...) =
    [detect_singularities(view(phase, :, :, :, i); keyargs...) for i in 1:size(phase, 4)]

residues(phase::LowDim{<:Real}; mask=nothing) = residues(to3d(phase); mask=to3d(mask))
detect_singularities(phase::LowDim{<:Real}; mask=nothing) = detect_singularities(to3d(phase); mask=to3d(mask))

## loop extraction
struct UnionFind
    parent::Vector{Int}
end
UnionFind(n::Integer) = UnionFind(collect(1:n))
function findroot(uf::UnionFind, x)
    while uf.parent[x] != x
        uf.parent[x] = uf.parent[uf.parent[x]]
        x = uf.parent[x]
    end
    return x
end
function unite!(uf::UnionFind, a, b)
    ra, rb = findroot(uf, a), findroot(uf, b)
    if ra != rb
        uf.parent[rb] = ra
    end
end

# position of `value` in the sorted vector `sorted`, 0 if not present
function sorted_lookup(sorted, value)
    i = searchsortedfirst(sorted, value)
    return (i ≤ length(sorted) && sorted[i] == value) ? i : 0
end

iscompletecube(J, sz) = all(d -> 1 ≤ J[d] < sz[d], 1:3)

# divergence of the residue field over the cube with the corner J
function cube_divergence(res, J)
    div = 0
    for d in 1:3
        e = unitindex(d)
        if checkbounds(Bool, res, d, J)
            div -= res[d, J]
        end
        if checkbounds(Bool, res, d, J + e)
            div += res[d, J+e]
        end
    end
    return div
end

# Two residue plaquettes belong to the same vortex line if they are faces of a
# common voxel cube. The divergence free residue field guarantees that the
# lines are closed as long as they don't leave the volume or the mask.
function find_loops(res::AbstractArray{Int8,4})
    sz = size(res)[2:4]
    lin = LinearIndices(res)
    linvox = LinearIndices(sz)
    idx = findall(!iszero, res)
    plaquettes = [lin[K] for K in idx] # ascending, `findall` returns linear order
    uf = UnionFind(length(idx))
    for (k, K) in enumerate(idx)
        c = K[1]
        I = CartesianIndex(K[2], K[3], K[4])
        for J in (I, I - unitindex(c)) # the two cubes sharing this plaquette
            all(>(0), Tuple(J)) || continue
            for d in 1:3, F in (J, J + unitindex(d)) # the 6 faces of the cube
                checkbounds(Bool, res, d, F) || continue
                iszero(res[d, F]) && continue # cheap test before the lookup
                q = sorted_lookup(plaquettes, lin[d, F])
                if q != 0
                    unite!(uf, k, q)
                end
            end
        end
    end

    groups = Dict{Int,Vector{Int}}()
    for k in eachindex(idx)
        push!(get!(groups, findroot(uf, k), Int[]), k)
    end

    loops = SingularityLoop[]
    for ks in values(groups)
        voxels = Set{Int}()
        closed = true
        for k in ks
            K = idx[k]
            c = K[1]
            I = CartesianIndex(K[2], K[3], K[4])
            a, b = PLAQUETTE_DIMS[c]
            ea, eb = unitindex(a), unitindex(b)
            for V in (I, I + ea, I + eb, I + ea + eb)
                push!(voxels, linvox[V])
            end
            for J in (I, I - unitindex(c))
                if !iscompletecube(J, sz) || cube_divergence(res, J) != 0
                    closed = false
                end
            end
        end
        push!(loops, SingularityLoop(sort!([lin[idx[k]] for k in ks]), sort!(collect(voxels)), closed))
    end
    sort!(loops; by=l -> (-loop_length(l), l.plaquettes[1])) # deterministic order
    return loops
end

## masks
"""
    singularity_mask(s::Singularities; min_length=1, max_length=Inf, closed_only=false)

Voxel mask of all voxels that are touched by a vortex line. `min_length` and
`max_length` filter on the number of residue plaquettes of a line. Loops with a
length of 4 are the smallest possible ones and are usually harmless noise.
"""
function singularity_mask(s::Singularities; min_length=1, max_length=Inf, closed_only=false)
    m = falses(s.size)
    for l in s.loops
        min_length ≤ loop_length(l) ≤ max_length || continue
        closed_only && !l.closed && continue
        m[l.voxels] .= true
    end
    return m
end

"""
    branchcuts(unwrapped::AbstractArray{<:Real,3}; threshold=π, mask=nothing)

Faces (connections between neighboring voxels) at which the unwrapped phase
jumps by more than `threshold`. These are the remaining branch cuts of the
congruent unwrapping, the surfaces that are bounded by the vortex lines.

`branchcuts(unwrapped)[d,i,j,k]` refers to the face between `(i,j,k)` and its
neighbor in the dimension `d`. Note that large true field jumps (air-tissue
interfaces) also show up here.
"""
function branchcuts(unwrapped::AbstractArray{<:Real,3}; threshold=Float64(π), mask=nothing)
    sz = size(unwrapped)
    mask = tomask(mask)
    cuts = falses(3, sz...)
    for d in 1:3
        e = unitindex(d)
        ranges = ntuple(t -> 1:(sz[t] - (t == d)), 3)
        for I in CartesianIndices(ranges)
            if !isnothing(mask) && !(mask[I] && mask[I+e])
                continue
            end
            jump = unwrapped[I+e] - unwrapped[I]
            cuts[d, I] = isfinite(jump) && abs(jump) > threshold
        end
    end
    return cuts
end

branchcuts(unwrapped::LowDim{<:Real}; mask=nothing, keyargs...) =
    branchcuts(to3d(unwrapped); mask=to3d(mask), keyargs...)

"""
    facemask_to_voxelmask(faces)

Reduce a face mask of size `(3, size...)` (as returned by [`branchcuts`](@ref))
to the mask of the voxels that are connected by these faces.
"""
function facemask_to_voxelmask(faces::AbstractArray{Bool,4})
    sz = size(faces)[2:4]
    m = falses(sz)
    for d in 1:3
        e = unitindex(d)
        ranges = ntuple(t -> 1:(sz[t] - (t == d)), 3)
        for I in CartesianIndices(ranges)
            if faces[d, I]
                m[I] = true
                m[I+e] = true
            end
        end
    end
    return m
end

"""
    congruence_error(corrected, original)

Deviation from congruence (`Δ ∈ 2πℤ`) between two phase images in radians.
Zero everywhere means that the corrected phase is still a valid unwrapping of
the input. Useful as a protection metric for [`fix_singularities!`](@ref).
"""
congruence_error(corrected, original) = abs.(rem2pi.(corrected .- original, RoundNearest))

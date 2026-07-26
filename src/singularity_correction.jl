## Correction of phase singularities
#
# ROMEO is congruent (the output differs from the input by n2π in every voxel),
# which is exactly the reason why a residue cannot disappear: the minimum
# spanning tree only moves the branch cut to the worst weights, it cannot
# remove it. Around a vein the vortex line encircles the vessel, so the cut
# surface has to leave the vessel somewhere into good tissue.
#
# The compromise implemented here is: congruent everywhere, except in small,
# explicitly marked patches. Every mode uses the same local weighted
# least-squares solve
#
#   min_φ Σ_ij W_ij (φ_i - φ_j - D_ij)²    with Dirichlet boundary conditions
#
# on a box around the vortex line, with the ROMEO result as boundary values.
# The least-squares solution keeps only the rotation free part of the field, so
# the branch cut is gone; outside of the patch nothing changes at all. A
# residue can only survive if the solution itself still changes by more than π
# between two voxels. The modes only differ in how the edge confidence w and
# the target gradient Δ are computed:
#
#   :lsq     w = mag², Δ from the measured wrapped phase differences
#   :smooth  w and Δ from the complex signal smoothed with `sigma`
#   :inpaint w = 0, Δ = 0 -> harmonic inpainting
#   :mask    no correction, only detection
#   :cascade small loops with :smooth, larger loops with :lsq (recommended)
#
# Writing W = w + ε and D = wΔ/(w+ε) turns the objective into
# Σ w(φ_i-φ_j-Δ)² + ε Σ (φ_i-φ_j)², so w -> 0 degenerates to inpainting
# without a second code path, and voxels with good signal (large w) keep their
# measured phase gradients and therefore stay congruent.

const SINGULARITY_MODES = (:off, :mask, :smooth, :lsq, :inpaint, :cascade)

"""
    fix_singularities!(phase; mode=:cascade, mag=nothing, mask=nothing, keyargs...)
    fix_singularities(phase; keyargs...) -> (phase, info)

Detect phase singularities in the unwrapped `phase` and remove them locally.
EXPERIMENTAL.

The unwrapped phase is only modified inside small patches around the detected
vortex lines, everywhere else it stays bit-identical (and therefore congruent
to the input phase). Returns a `NamedTuple` with the detection and correction
outputs, `phase` is modified in place.

### Modes

- `:off`: do nothing
- `:mask`: only detect, don't touch the phase. The mask can be used as a data
    term weight for a downstream QSM inversion.
- `:cascade`: (recommended) `:smooth` for small loops, `:lsq` for larger ones
- `:smooth`: complex smoothing with `sigma` annihilates ± pairs inside the
    kernel width. Only suitable for small loops, the winding number of a large
    loop is topologically protected and would only be moved.
- `:lsq`: magnitude weighted least-squares solve in the patch
- `:inpaint`: harmonic inpainting in the patch (discards the intravascular
    phase, use as a control arm)

### Optional keyword arguments

- `mag`: magnitude, used for the edge weights and for the gate
- `mask`: only detect and correct inside the mask
- `sigma=1.0`: kernel width in voxels for the complex smoothing
- `small_loop_length=8`: loops up to this number of residue plaquettes are
    treated as small in the `:cascade` mode
- `min_loop_length=1`: shorter loops are not corrected (4 is the smallest
    possible loop and is usually noise)
- `max_loop_length=1000`: longer loops are only marked, never touched. This is
    the fallback that protects real pathology (microbleeds, implants).
- `max_patch_voxels=20000`: larger patches are only marked
- `dilate=2`: dilation of the vortex line and its branch cut surface, defines
    the corrected region
- `pad=2`: width of the Dirichlet boundary ring around the patch
- `mag_threshold=0.2`: magnitude relative to its 95th percentile below which a
    voxel counts as a signal void
- `min_void_voxels=1`: number of signal void voxels a vortex line has to touch
    to be corrected at all. A singularity is a zero of the complex signal, a
    line without a void is a different phenomenon (e.g. an aliased field
    gradient in good tissue) and is only marked.
- `require_closed=false`: if `true`, only lines that close inside the volume
    are corrected. Lines that leave the volume or the mask are consistent as
    well, as long as the patch contains everything that stays inside.
- `epsilon=1e-3`: regularization of the least-squares problem
- `cut_threshold=π`: phase jump above which a face counts as a branch cut

### Output fields

- `singularities`: voxel mask of all detected vortex lines
- `cuts`: voxel mask of the remaining branch cuts (jumps > `cut_threshold`)
- `changed`: voxel mask of the actually corrected voxels
- `delta`: the correction `φ_corrected - φ_ROMEO` (empty for `:mask`/`:off`)
- `loops`, `nloops`, `nresidues`, `nchanged`
- `nfixed`/`nskipped`: number of vortex lines that were corrected / only marked

# Examples
```julia-repl
julia> unwrapped = unwrap(phase; mag, TEs);
julia> info = fix_singularities!(unwrapped; mag, mode=:cascade);
julia> info.nfixed, info.nskipped
julia> savenii(info.singularities, "singularities.nii"; header=header(phase));
```

See also [`detect_singularities`](@ref), [`branchcuts`](@ref)
"""
fix_singularities, fix_singularities!

function fix_singularities(phase; keyargs...)
    phase = copy(phase)
    return phase, fix_singularities!(phase; keyargs...)
end

function fix_singularities!(phase::AbstractArray{T,3};
        mode=:cascade, mag=nothing, mask=nothing,
        sigma=1.0, small_loop_length=8, min_loop_length=1, max_loop_length=1000,
        max_patch_voxels=20_000, dilate=2, pad=2, mag_threshold=0.2,
        min_void_voxels=1, require_closed=false, epsilon=1e-3,
        cut_threshold=Float64(π), keyargs...) where T
    mode = Symbol(mode)
    mode in SINGULARITY_MODES || throw(ArgumentError("singularity mode '$mode' is undefined! Use one of $(join(SINGULARITY_MODES, ", "))"))
    sz = size(phase)
    for (name, arr) in (("mag", mag), ("mask", mask))
        if !isnothing(arr) && size(arr) != sz
            throw(DimensionMismatch("size of $name is $(size(arr)), but the phase is $sz"))
        end
    end
    correcting = mode ∉ (:off, :mask)
    if mode == :off
        return (; mode, singularities=falses(sz), cuts=falses(sz), changed=falses(sz),
            delta=zeros(T, 0, 0, 0), loops=SingularityLoop[], nloops=0, nfixed=0,
            nskipped=0, nresidues=0, nchanged=0)
    end

    mask = tomask(mask)
    sing = detect_singularities(phase; mask)
    cuts = branchcuts(phase; threshold=cut_threshold, mask)
    changed = falses(sz)
    delta = zeros(T, correcting ? sz : (0, 0, 0))
    handled = falses(length(sing.loops))
    cutmask = facemask_to_voxelmask(cuts)

    if correcting
        maxmag = magnitude_scale(mag)
        gate = singularity_gate(mag, mask, mag_threshold, maxmag, sz)
        allowed = isnothing(mask) ? trues(sz) : mask
        plaq_index = plaquette_index(sing.loops)
        cuts_components = cut_components(cutmask, allowed, sz)
        on_loop = singularity_mask(sing)
        for (i, loop) in enumerate(sing.loops)
            handled[i] && continue
            is_correctable(loop, gate; min_loop_length, max_loop_length, min_void_voxels, require_closed) || continue
            patch = build_patch(loop, sing, cuts_components, on_loop, plaq_index, sz; dilate, pad, max_patch_voxels)
            isnothing(patch) && continue
            box_lo, box_hi, interior, ids = patch
            # A line that is protected by the length limit must not be modified
            # just because it shares a patch with a correctable one.
            any(id -> loop_length(sing.loops[id]) > max_loop_length, ids) && continue
            lmode = if mode == :cascade
                loop_length(loop) ≤ small_loop_length ? :smooth : :lsq
            else
                mode
            end
            fix_patch!(phase, delta, changed, box_lo, box_hi, interior, mag, mask, lmode;
                sigma, epsilon, maxmag)
            for id in ids
                handled[id] = true
            end
        end
    end

    return (; mode,
        singularities=singularity_mask(sing),
        cuts=cutmask,
        changed, delta,
        loops=sing.loops,
        nloops=length(sing.loops),
        nfixed=count(handled),
        nskipped=length(sing.loops) - count(handled),
        nresidues=sum(abs, sing.residues),
        nchanged=count(changed))
end

fix_singularities!(phase::LowDim{<:Real}; mag=nothing, mask=nothing, keyargs...) =
    fix_singularities!(to3d(phase); mag=to3d(mag), mask=to3d(mask), keyargs...)

function fix_singularities!(phase::AbstractArray{T,4}; mag=nothing, keyargs...) where T
    infos = map(1:size(phase, 4)) do i
        m = isnothing(mag) ? nothing : view(mag, :, :, :, min(i, size(mag, 4)))
        fix_singularities!(view(phase, :, :, :, i); mag=m, keyargs...)
    end
    stack4(f) = cat((getfield(info, f) for info in infos)...; dims=4)
    return (; mode=infos[1].mode,
        singularities=stack4(:singularities),
        cuts=stack4(:cuts),
        changed=stack4(:changed),
        delta=stack4(:delta),
        loops=[info.loops for info in infos],
        nloops=sum(info -> info.nloops, infos),
        nfixed=sum(info -> info.nfixed, infos),
        nskipped=sum(info -> info.nskipped, infos),
        nresidues=sum(info -> info.nresidues, infos),
        nchanged=sum(info -> info.nchanged, infos),
        echoes=infos)
end

## gating: where are we allowed to touch the phase
magnitude_scale(mag) = isnothing(mag) ? 1.0 : max(eps(), Float64(quantile(filter(isfinite, mag), 0.95)))

# The physical requirement is "don't change the phase where there is signal".
# A singularity sits where the complex signal vanishes, so a vortex line has to
# touch a signal void to be the case we are after. The magnitude weighting of
# the least-squares solve enforces the same requirement a second time inside
# the patch.
function singularity_gate(mag, mask, mag_threshold, maxmag, sz)
    gate = if isnothing(mag)
        trues(sz)
    else
        BitArray(.!isfinite.(mag) .| (mag .< mag_threshold * maxmag))
    end
    if !isnothing(mask)
        gate .&= mask
    end
    return gate
end

# A vortex line that leaves the volume or the mask is not closed, but it can
# still be corrected: nothing crosses the patch boundary there. Only the part
# of a line that stays inside must be enclosed by the patch, which is checked
# in `build_patch`.
function is_correctable(loop, gate; min_loop_length, max_loop_length, min_void_voxels, require_closed)
    require_closed && !loop.closed && return false
    min_loop_length ≤ loop_length(loop) ≤ max_loop_length || return false
    return count(v -> gate[v], loop.voxels) ≥ min_void_voxels
end

## patch construction
# plaquette linear index -> loop number, as a sorted vector pair
function plaquette_index(loops)
    keys = Int[]
    ids = Int[]
    for (i, l) in enumerate(loops), p in l.plaquettes
        push!(keys, p)
        push!(ids, i)
    end
    perm = sortperm(keys)
    return keys[perm], ids[perm]
end

function build_patch(loop, sing, cuts_components, on_loop, plaq_index, sz; dilate, pad, max_patch_voxels)
    seed = cut_seed(loop, cuts_components, sz, max_patch_voxels)
    isnothing(seed) && return nothing
    interior = dilate_voxels(seed, dilate, sz, max_patch_voxels)
    isnothing(interior) && return nothing

    lo, hi = voxel_bounds(interior, sz)
    box_lo = CartesianIndex(ntuple(d -> max(1, lo[d] - pad), 3))
    box_hi = CartesianIndex(ntuple(d -> min(sz[d], hi[d] + pad), 3))

    # The net charge crossing the patch boundary has to be zero, otherwise the
    # boundary condition is inconsistent and the error is smeared over the
    # patch. Every vortex line touching the patch must therefore be enclosed.
    ids = touching_loops(interior, sing, on_loop, plaq_index, sz)
    for id in ids
        all(v -> v in interior, sing.loops[id].voxels) || return nothing
    end
    return box_lo, box_hi, interior, ids
end

# The branch cut surface is bounded by the vortex line, but the minimum
# spanning tree places it wherever the weights are worst, not necessarily next
# to the line. The patch has to cover the whole surface, otherwise the 2π jump
# survives outside of it. The connected components of the cut voxels are
# labeled once for the whole volume, a patch then consists of the vortex line
# plus the components it touches.
function cut_components(cutmask, allowed, sz)
    labels = zeros(Int32, sz)
    voxels = Vector{Vector{Int}}()
    car = CartesianIndices(sz)
    lin = LinearIndices(sz)
    queue = Int[]
    for v in eachindex(labels)
        (cutmask[v] && allowed[v] && iszero(labels[v])) || continue
        push!(voxels, Int[])
        label = Int32(length(voxels))
        labels[v] = label
        empty!(queue)
        push!(queue, v)
        while !isempty(queue)
            w = pop!(queue)
            push!(voxels[label], w)
            I = car[w]
            for d in 1:3, s in (-1, 1)
                J = I + s * unitindex(d)
                checkbounds(Bool, labels, J) || continue
                (cutmask[J] && allowed[J] && iszero(labels[J])) || continue
                labels[J] = label
                push!(queue, lin[J])
            end
        end
    end
    return labels, voxels
end

function cut_seed(loop, (labels, comp_voxels), sz, limit)
    car = CartesianIndices(sz)
    components = Set{Int32}()
    for v in loop.voxels
        iszero(labels[v]) || push!(components, labels[v])
        I = car[v]
        for d in 1:3, s in (-1, 1)
            J = I + s * unitindex(d)
            checkbounds(Bool, labels, J) || continue
            iszero(labels[J]) || push!(components, labels[J])
        end
    end
    total = length(loop.voxels) + sum(c -> length(comp_voxels[c]), components; init=0)
    total > limit && return nothing
    seed = Set{Int}(loop.voxels)
    for c in components, v in comp_voxels[c]
        push!(seed, v)
    end
    return seed
end

function dilate_voxels(seed, r, sz, limit)
    out = Set{Int}()
    car = CartesianIndices(sz)
    lin = LinearIndices(sz)
    for v in seed
        I = car[v]
        for J in CartesianIndices(ntuple(d -> max(1, I[d] - r):min(sz[d], I[d] + r), 3))
            push!(out, lin[J])
        end
        length(out) > limit && return nothing
    end
    return out
end

function voxel_bounds(voxels, sz)
    car = CartesianIndices(sz)
    lo = [typemax(Int), typemax(Int), typemax(Int)]
    hi = [0, 0, 0]
    for v in voxels
        I = car[v]
        for d in 1:3
            lo[d] = min(lo[d], I[d])
            hi[d] = max(hi[d], I[d])
        end
    end
    return CartesianIndex(lo...), CartesianIndex(hi...)
end

function touching_loops(interior, sing, on_loop, (plaq_keys, plaq_ids), sz)
    ids = Set{Int}()
    car = CartesianIndices(sz)
    lin = LinearIndices(sing.residues)
    for v in interior
        on_loop[v] || continue # only these can be a corner of a residue plaquette
        I = car[v]
        for c in 1:3
            a, b = PLAQUETTE_DIMS[c]
            ea, eb = unitindex(a), unitindex(b)
            for P in (I, I - ea, I - eb, I - ea - eb)
                all(>(0), Tuple(P)) || continue
                (P[a] < sz[a] && P[b] < sz[b]) || continue
                sing.residues[c, P] == 0 && continue
                j = sorted_lookup(plaq_keys, lin[c, P])
                j != 0 && push!(ids, plaq_ids[j])
            end
        end
    end
    return ids
end

## local weighted least-squares solve
function fix_patch!(phase, delta, changed, box_lo, box_hi, interior, mag, mask, lmode; sigma, epsilon, maxmag)
    bsz = Tuple(box_hi - box_lo) .+ 1
    lin = LinearIndices(size(phase))
    φ = zeros(Float64, bsz)
    free = falses(bsz)
    valid = falses(bsz)
    for L in CartesianIndices(bsz)
        I = box_lo + L - oneunit(L)
        p = Float64(phase[I])
        φ[L] = isfinite(p) ? p : 0.0
        valid[L] = isfinite(p) && (isnothing(mask) || mask[I])
        free[L] = valid[L] && (lin[I] in interior)
    end
    any(free) || return

    W = zeros(Float64, 3, bsz...)
    D = zeros(Float64, 3, bsz...)
    edge_terms!(W, D, phase, mag, box_lo, valid, lmode, sigma, maxmag, epsilon)
    solve_lsq!(φ, free, W, D)

    for L in CartesianIndices(bsz)
        (free[L] && isfinite(φ[L])) || continue
        I = box_lo + L - oneunit(L)
        value = eltype(phase)(φ[L])
        value == phase[I] && continue
        delta[I] = value - phase[I]
        phase[I] = value
        changed[I] = true
    end
    return
end

function edge_terms!(W, D, phase, mag, box_lo, valid, lmode, sigma, maxmag, epsilon)
    bsz = size(valid)
    box_hi = box_lo + CartesianIndex(bsz) - oneunit(box_lo)
    S, offset = if lmode == :smooth
        s, slo = smoothed_signal(phase, mag, box_lo, box_hi, sigma)
        s, box_lo - slo
    else
        nothing, CartesianIndex(0, 0, 0)
    end
    for d in 1:3
        e = unitindex(d)
        for L in CartesianIndices(ntuple(t -> 1:(bsz[t] - (t == d)), 3))
            M = L + e
            (valid[L] && valid[M]) || continue
            I = box_lo + L - oneunit(L)
            J = I + e
            w, Δ = if lmode == :inpaint
                0.0, 0.0
            elseif lmode == :smooth
                sI, sJ = S[L+offset], S[M+offset]
                z = sJ * conj(sI)
                abs(sI) * abs(sJ) / maxmag^2, iszero(z) ? 0.0 : angle(z)
            else
                m = isnothing(mag) ? 1.0 : Float64(mag[I]) * Float64(mag[J]) / maxmag^2
                m, rem2pi(Float64(phase[J]) - Float64(phase[I]), RoundNearest)
            end
            isfinite(w) || (w = 0.0)
            isfinite(Δ) || (Δ = 0.0)
            w = clamp(w, 0.0, 1.0)
            W[d, L] = w + epsilon
            D[d, L] = w * Δ / (w + epsilon)
        end
    end
    return
end

# Conjugate gradient on the weighted graph Laplacian with Dirichlet boundary
# conditions. `free` marks the unknowns, all other voxels keep their value.
function solve_lsq!(φ, free, W, D; maxiter=2000, tol=1e-8)
    bsz = size(φ)
    Wsum = zeros(Float64, bsz)
    for d in 1:3
        e = unitindex(d)
        for I in CartesianIndices(ntuple(t -> 1:(bsz[t] - (t == d)), 3))
            w = W[d, I]
            Wsum[I] += w
            Wsum[I+e] += w
        end
    end
    free = free .& (Wsum .> 0)
    r = zeros(Float64, bsz)
    q = zeros(Float64, bsz)
    lsq_residual!(r, φ, free, W, D, Wsum)
    p = copy(r)
    rs = dotp(r, r)
    rs0 = rs
    iszero(rs0) && return 0.0
    for _ in 1:maxiter
        lsq_apply!(q, p, free, W, Wsum)
        pq = dotp(p, q)
        pq ≤ 0 && break
        α = rs / pq
        @. φ += α * p # p is zero outside of the free voxels
        @. r -= α * q
        rsnew = dotp(r, r)
        if rsnew ≤ tol^2 * rs0
            rs = rsnew
            break
        end
        @. p = r + (rsnew / rs) * p
        rs = rsnew
    end
    return sqrt(rs / rs0)
end

function lsq_residual!(r, φ, free, W, D, Wsum)
    bsz = size(φ)
    fill!(r, 0)
    for d in 1:3
        e = unitindex(d)
        for I in CartesianIndices(ntuple(t -> 1:(bsz[t] - (t == d)), 3))
            w = W[d, I]
            iszero(w) && continue
            J = I + e
            r[I] += w * (φ[J] - D[d, I])
            r[J] += w * (φ[I] + D[d, I])
        end
    end
    @. r = (r - φ * Wsum) * free
    return r
end

function lsq_apply!(q, p, free, W, Wsum)
    bsz = size(p)
    fill!(q, 0)
    for d in 1:3
        e = unitindex(d)
        for I in CartesianIndices(ntuple(t -> 1:(bsz[t] - (t == d)), 3))
            w = W[d, I]
            iszero(w) && continue
            J = I + e
            q[I] -= w * p[J]
            q[J] -= w * p[I]
        end
    end
    @. q = (q + p * Wsum) * free
    return q
end

function dotp(a, b)
    s = 0.0
    @inbounds @simd for i in eachindex(a, b)
        s += a[i] * b[i]
    end
    return s
end

## complex smoothing
function smoothed_signal(phase, mag, lo, hi, sigma)
    sz = size(phase)
    ext = max(1, ceil(Int, 3sigma))
    elo = CartesianIndex(ntuple(d -> max(1, lo[d] - ext), 3))
    ehi = CartesianIndex(ntuple(d -> min(sz[d], hi[d] + ext), 3))
    esz = Tuple(ehi - elo) .+ 1
    S = zeros(ComplexF64, esz)
    for L in CartesianIndices(esz)
        I = elo + L - oneunit(L)
        m = isnothing(mag) ? 1.0 : Float64(mag[I])
        p = Float64(phase[I])
        if isfinite(m) && isfinite(p)
            S[L] = m * cis(p)
        end
    end
    smooth_complex!(S, sigma)
    return S, elo
end

function gaussiankernel(sigma; truncate=3)
    r = max(1, ceil(Int, truncate * sigma))
    k = [exp(-0.5(i / sigma)^2) for i in -r:r]
    return k ./ sum(k)
end

function smooth_complex!(S, sigma)
    k = gaussiankernel(sigma)
    r = (length(k) - 1) ÷ 2
    tmp = similar(S)
    for d in 1:3
        size(S, d) == 1 && continue
        convolve_dim!(tmp, S, k, r, d)
        copyto!(S, tmp)
    end
    return S
end

function convolve_dim!(out, A, k, r, d)
    sz = size(A)
    for I in CartesianIndices(A)
        acc = zero(eltype(A))
        for (n, w) in enumerate(k)
            j = clamp(I[d] + n - r - 1, 1, sz[d])
            acc += w * A[CartesianIndex(ntuple(t -> t == d ? j : I[t], 3))]
        end
        out[I] = acc
    end
    return out
end

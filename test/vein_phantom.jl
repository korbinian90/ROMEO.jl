## Numerical vein phantom
#
# A cylindrical vessel with the susceptibility difference `dchi` at the angle
# `theta` to B0 is placed in a box. The analytic field of an infinite cylinder
# is evaluated on a `sub`^3 subgrid, the complex signal (with T2* decay) is
# generated there and averaged over the subvoxels of each voxel.
#
# Nothing is put in by hand: the partial voluming at the vessel wall and the
# intravoxel dephasing in the steep dipole field make the complex signal cross
# zero somewhere, which is exactly the mechanism that creates phase
# singularities in high resolution in-vivo data. The phantom therefore delivers
# singularities together with the true field map.

using Random

const γ_over_2π = 42.577478e6 # Hz/T

"""
    vein_phantom(; keyargs...)

Returns `(; phase, mag, complex, field, vein)` for a single straight vessel.
`field` is the true frequency offset in Hz, `vein` the true intravascular mask
(both on the voxel grid), `phase` is wrapped to [-π;π].
"""
function vein_phantom(;
        sz=(40, 40, 24),
        voxelsize=0.3,          # mm, isotropic
        radius=0.45,            # mm, vessel radius
        theta=deg2rad(70),      # angle between vessel and B0
        dchi=0.4e-6,            # susceptibility difference vein - tissue
        B0=7.0,                 # T
        TE=25e-3,               # s
        T2s_tissue=25e-3,
        T2s_blood=8e-3,
        rho_tissue=1.0,
        rho_blood=1.0,
        sub=3,                  # subvoxels per dimension
        noise=0.0,
        seed=0,
        )
    rng = MersenneTwister(seed)
    axis = (sin(theta), 0.0, cos(theta))            # vessel direction
    # projection of B0 (z) into the plane perpendicular to the vessel
    zperp = (0.0, 0.0, 1.0) .- cos(theta) .* axis
    nzperp = sqrt(sum(abs2, zperp))
    zperp = nzperp > 0 ? zperp ./ nzperp : (1.0, 0.0, 0.0)

    center = (sz .+ 1) ./ 2
    S = zeros(ComplexF64, sz)
    field = zeros(Float64, sz)
    vein = falses(sz)
    step = voxelsize / sub
    offsets = ((1:sub) .- (sub + 1) / 2) .* step

    for I in CartesianIndices(sz)
        acc = zero(ComplexF64)
        facc = 0.0
        inside_count = 0
        for dx in offsets, dy in offsets, dz in offsets
            r = ((I[1] - center[1]) * voxelsize + dx,
                 (I[2] - center[2]) * voxelsize + dy,
                 (I[3] - center[3]) * voxelsize + dz)
            f, inside = cylinder_frequency(r, axis, zperp, radius, dchi, theta, B0)
            inside_count += inside
            facc += f
            m = inside ? rho_blood * exp(-TE / T2s_blood) : rho_tissue * exp(-TE / T2s_tissue)
            acc += m * cis(2π * f * TE)
        end
        n = sub^3
        S[I] = acc / n
        field[I] = facc / n
        vein[I] = 2inside_count > n
    end
    if noise > 0
        S .+= noise .* complex.(randn(rng, sz), randn(rng, sz))
    end
    return (; phase=angle.(S), mag=abs.(S), complex=S, field, vein)
end

# analytic field of an infinite cylinder (Lorentz sphere corrected)
function cylinder_frequency(r, axis, zperp, radius, dchi, theta, B0)
    along = sum(r .* axis)
    perp = r .- along .* axis
    ρ = sqrt(sum(abs2, perp))
    if ρ ≤ radius
        return γ_over_2π * B0 * dchi / 6 * (3cos(theta)^2 - 1), true
    end
    cosφ = sum(perp .* zperp) / ρ
    cos2φ = 2cosφ^2 - 1
    return γ_over_2π * B0 * dchi / 2 * (radius / ρ)^2 * sin(theta)^2 * cos2φ, false
end

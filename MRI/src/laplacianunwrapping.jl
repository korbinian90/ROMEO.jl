    using FFTW

    laplacianunwrap(ϕ) = laplacianunwrap!(copy(ϕ))
    function laplacianunwrap!(ϕ)
        FFTW.set_num_threads(Threads.nthreads())
        ϕ .+= 2π .* k(ϕ)
    end

    # Schofield and Zhu 2003, https://doi.org/10.1364/OL.28.001194
    k(ϕw) = 1 / 2π .* ∇⁻²(∇²_nw(ϕw) - ∇²(ϕw))  # (1)

    ∇²(x) = -(2π)^ndims(x) / length(x) .* idct(pqterm(size(x)) .* dct(x))  # (2)

    ∇⁻²(x) = -length(x) / (2π)^ndims(x) .* idct(dct(x) ./ pqterm(size(x)))  # (3)

    ∇²_nw(ϕw) = cos.(ϕw) .* ∇²(sin.(ϕw)) .- sin.(ϕw) .* ∇²(cos.(ϕw))  # (in text)


    pqterm(sz::NTuple{1}) = (1:sz[1]).^2  # 1D case
    pqterm(sz::NTuple{2}) = [p^2 + q^2 for p in 1:sz[1], q in 1:sz[2]]  # 2D case
    pqterm(sz::NTuple{3}) = [p^2 + q^2 + t^2 for p in 1:sz[1], q in 1:sz[2], t in 1:sz[3]]  # 3D case
    pqterm(sz::NTuple{4}) = [p^2 + q^2 + t^2 + r^2 for p in 1:sz[1], q in 1:sz[2], t in 1:sz[3], r in 1:sz[4]]  # 4D case

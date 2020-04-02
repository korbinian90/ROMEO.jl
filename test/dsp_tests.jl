@testset "Unwrap 1D" begin
    @test unwrap([0.1, 0.2, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 - 2pi, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 - 2pi, 0.3 - 2pi, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1 + 2pi, 0.2, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]

    test_v = [0.1, 0.2, 0.3 + 2pi, 0.4]
    res_v = unwrap(test_v)
    @test test_v ≈ [0.1, 0.2, 0.3 + 2pi, 0.4]
    res_v .= 0
    unwrap!(test_v)
    @test test_v ≈ [0.1, 0.2, 0.3, 0.4]

    # test unwrapping within multi-dimensional array
    wrapped = [0.1, 0.2 + 2pi, 0.3, 0.4]
    unwrapped = [0.1, 0.2, 0.3, 0.4]
    wrapped = hcat(wrapped, wrapped)
    unwrapped = hcat(unwrapped, unwrapped)

    # test generically typed unwrapping
    types = (Float32, Float64)
    for T in types
        A_unwrapped = collect(range(0, stop=4convert(T, π), length=10))
        A_wrapped = A_unwrapped .% (2convert(T, π))

        test(I, J) = I ≈ J || I .+ 2π ≈ J || I .+ 4π ≈ J

        @test test(unwrap(A_wrapped), A_unwrapped)
        unwrap!(A_wrapped)
        @test test(A_wrapped, A_unwrapped)
    end
end

# tests for multi-dimensional unwrapping
@testset "Unwrap 2D" begin
    types = (Float32, Float64)
    for T in types
        v_unwrapped = collect(range(0, stop=4convert(T, π), length=7))
        A_unwrapped = v_unwrapped .+ v_unwrapped'
        A_wrapped = A_unwrapped .% (2convert(T, π))

        test_unwrapped = unwrap(A_wrapped)
        d = first(A_unwrapped) - first(test_unwrapped)
        @test (test_unwrapped .+ d) ≈ A_unwrapped
        unwrap!(A_wrapped)
        d = first(A_unwrapped) - first(A_wrapped)
        @test (A_wrapped .+ d) ≈ A_unwrapped

        v_unwrapped_range = collect(range(0, stop=4, length=7))
        A_unwrapped_range = v_unwrapped_range .+ v_unwrapped_range'
        test_range = convert(T, 2)
        A_wrapped_range = A_unwrapped_range .% test_range

    end
end

@testset "Unwrap 3D" begin
    types = (Float32, Float64)
    f(x, y, z) = 0.1x^2 - 2y + 2z
    f_wraparound2(x, y, z) = 5*sin(x) + 2*cos(y) + z
    f_wraparound3(x, y, z) = 5*sin(x) + 2*cos(y) - 4*cos(z)
    for T in types
        grid = range(zero(T), stop=2convert(T, π), length=11)
        f_uw = f.(grid, grid', reshape(grid, 1, 1, :))
        f_wr = f_uw .% (2convert(T, π))
        uw_test = unwrap(f_wr)
        offset = first(f_uw) - first(uw_test)
        @test (uw_test.+offset) ≈ f_uw rtol=eps(T) #oop, nowrap
        # test in-place version
        unwrap!(f_wr)
        offset = first(f_uw) - first(f_wr)
        @test (f_wr.+offset) ≈ f_uw rtol=eps(T) #ip, nowrap

        f_uw = f_wraparound2.(grid, grid', reshape(grid, 1, 1, :))
        f_wr = f_uw .% (2convert(T, π))
        uw_test = unwrap(f_wr)
        offset = first(f_uw) - first(uw_test)
        @test (uw_test.+offset) ≈ f_uw #oop, 2wrap
        # test in-place version
        unwrap!(f_wr)
        offset = first(f_uw) - first(f_wr)
        @test (f_wr.+offset) ≈ f_uw #ip, 2wrap

        f_uw = f_wraparound3.(grid, grid', reshape(grid, 1, 1, :))
        f_wr = f_uw .% (2convert(T, π))
        uw_test = unwrap(f_wr)
        offset = first(f_uw) - first(uw_test)
        @test (uw_test.+offset) ≈ f_uw #oop, 3wrap
        # test in-place version
        unwrap!(f_wr)
        offset = first(f_uw) - first(f_wr)
        @test (f_wr.+offset) ≈ f_uw #oop, 3wrap
    end
end

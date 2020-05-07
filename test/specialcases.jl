@testset "Special Cases" begin

phase = ones(3,3,3)
unwrap(phase)

point = ones(1,1,1)
@test_throws AssertionError unwrap(point)

## weight test
function weight_test(w, out)
    w = reshape(w, length(w), 1, 1)
    @test UInt8.(out) == ROMEO.calculateweights(w)[1,:,1,1]
end
@test UInt[8, 8, 8, 0] == ROMEO.calculateweights(reshape([0.1, 0.2 + 2pi, 0.3, 0.4], 4,1,1))[1,:,1,1]
weight_test([0.1, 0.2 + 2pi, 0.3, 0.4], [8, 8, 8, 0])
weight_test([0.1, 0.2 + 2pi, 0.3, NaN], [8, 255, 0, 0]) # 255 (worst valid value) at voxel bordering to NaN voxel

## NaN test
@test nan_test(unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test nan_test(unwrap([0.1, 0.2 + 2pi, 0.3, NaN]), [0.1, 0.2, 0.3, NaN])
@test nan_test(unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]; mag=[1, 1, 1, NaN]), [0.1, 0.2, 0.3, 0.4])

## test dimensions
function dim_test(a; keyargs...)
    ca = copy(a)
    ua = unwrap(a; keyargs...)
    @test a != ua
    @test size(a) == size(ua)
    @test a == ca
    unwrap!(a; keyargs...)
    @test a != ca
    @test a == ua
end

dim_test(10rand(50))
dim_test(10rand(50,20))
dim_test(10rand(15,20,10))
dim_test(10rand(7,9,10,3); TEs=ones(3))

end

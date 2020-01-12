v = -5:1.0:5
for a in [v, v .+ v', v .+ v' .+ reshape(v, 1,1,length(v))]
    wrapped = rem2pi.(a, RoundNearest)
    @test unwrapscore(a) < 0.1
    @test 0.1 < unwrapscore(a .+ 30) < 1
    @test unwrapscore(wrapped) > 5
end

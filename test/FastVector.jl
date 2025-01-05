using Test, PolynomialOptimization.FastVector

@testset "FastVector (behave as Vector)" begin
    # examples from Julia doc
    v = FastVec{Int}(undef, 3)
    v .= 1:3
    @test push!(v, 4, 5, 6) == [1, 2, 3, 4, 5, 6]

    v = FastVec{Any}(undef, 6)
    v .= 1:6
    @test insert!(v, 3, "here") == Any[1, 2, "here", 3, 4, 5, 6]

    v = FastVec{Int}(undef, 6)
    v .= 6:-1:1
    @test deleteat!(v, 2) == [6, 4, 3, 2, 1]

    v = FastVec{Int}(undef, 6) # we only support unit range
    v .= 6:-1:1
    @test deleteat!(v, 3:5) == [6, 5, 1]

    v = FastVec{Int}(undef, 6) # we only support unit range
    v .= 6:-1:1
    @test splice!(v, 5) == 2
    @test v == [6, 5, 4, 3, 1]
    @test splice!(v, 5, -1) == 1
    @test v == [6, 5, 4, 3, -1]
    @test splice!(v, 1, [-1, -2, -3]) == 6
    @test v == [-1, -2, -3, 5, 4, 3, -1]
    @test splice!(v, 4:3, 2) == Int[]
    @test v == [-1, -2, -3, 2, 5, 4, 3, -1]
    @test splice!(v, 2:3) == [-2, -3]
    @test v == [-1, 2, 5, 4, 3, -1]
    @test splice!(v, 2:4, [-2, -5, -4]) == [2, 5, 4]
    @test v == [-1, -2, -5, -4, 3, -1]
    @test splice!(v, 2:4, [1, 2, 3, 4, 5]) == [-2, -5, -4]
    @test v == [-1, 1, 2, 3, 4, 5, 3, -1]
    @test splice!(v, 2:4, [7]) == [1, 2, 3]
    @test v == [-1, 7, 4, 5, 3, -1]
    @test splice!(v, 4:0, [9, 11, 13]) == Int[]
    @test v == [-1, 7, 4, 9, 11, 13, 5, 3, -1]

    v = FastVec{Int}(undef, 1)
    v[1] = 1
    @test append!(v, [2, 3]) == [1, 2, 3]
    @test append!(v, [4, 5], [6]) == [1, 2, 3, 4, 5, 6]

    v = FastVec{Int}(undef, 1)
    v[1] = 3
    @test prepend!(v, [1, 2]) == [1, 2, 3]

    v = FastVec{Int}(undef, 1)
    v[1] = 6
    @test prepend!(v, [1, 2], [3, 4, 5]) == [1, 2, 3, 4, 5, 6]
end
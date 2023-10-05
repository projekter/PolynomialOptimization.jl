using Test

#optimize = all(x != "noopt" for x in ARGS);
optimize = true
@testset "PolynomialOptimization" begin
    @testset "no sparsity" begin
        include("./SparsityNone.jl")
    end
    @testset "correlative sparsity" begin
        include("./SparsityCorrelative.jl")
    end
    @testset "term sparsity" begin
        include("./SparsityTermBlock.jl")
    end
    @testset "term sparsity with chordal extension" begin
        include("./SparsityTermCliques.jl")
    end
    @testset "correlative and term sparsity" begin
        include("./SparsityCorrelativeTerm.jl")
    end
    @testset "tightening using Nie's method" begin
        include("./Tightening.jl")
    end
    @testset "multiplication by prefactor using Mai et al.'s method" begin
        include("./Noncompact.jl")
    end
    @testset "Packed matrices" begin
        include("./PackedMatrices.jl")
    end
    @testset "LANCELOT solver" begin
        include("./LocalSolver.jl")
    end
end
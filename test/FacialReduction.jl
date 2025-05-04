include("./shared.jl")
using LinearAlgebra
using PolynomialOptimization.FacialReduction: facial_reduction, reduction_methods

macro testfrs(fn)
    esc(quote
        for method in reduction_methods
            $fn
        end
    end)
end

@testset "Permenter, Parillo IV.C" begin
    DynamicPolynomials.@polyvar x[1:2]
    prob = poly_problem(x[1]^2*x[2]^4, zero=[1-x[1]^4*x[2]^2], soscert=true)
    @testfrs begin
        gr = Relaxation.groupings(Relaxation.FacialReduction(prob, 3))
        @test gr.obj[1] == IntMonomialVector([x[1]*x[2]^2])
    end
end

@testset "Permenter, Parillo V.A van der Pol Oscillator" begin
    DynamicPolynomials.@polyvar x[1:2]
    V = 15x[1]^2 - 10x[1]*x[2] + 10x[2]^2
    Vdot = -10x[1]^2 - 10x[1]^3*x[2] - 10x[2]^2 + 20x[1]^2*x[2]^2
    prob = poly_problem(V, factor_coercive=x[1]^2 + x[2]^2, zero=[Vdot])
    @testfrs begin
        gr = Relaxation.groupings(Relaxation.FacialReduction(prob, 4))
        @test gr.obj[1] == IntMonomialVector([x[2], x[2]^2, x[1], x[1]*x[2], x[1]*x[2]^2, x[1]^2, x[1]^2*x[2], x[1]^2*x[2]^2,
                                              x[1]^3*x[2]])
    end
end

# Unfortunately, we cannot test the second example, as the script to generate the matrices is no longer available.
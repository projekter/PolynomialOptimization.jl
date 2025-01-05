include("./shared.jl")

@testset "Example 3.1" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityCorrelativeTerm(poly_problem(1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] * x[2] + x[2] * x[3] + x[3]), 1)
    @test strRep(sp) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2]
  PSD block sizes:
    [2 => 1]
> Clique #2: x[2], x[3]
  PSD block sizes:
    [3 => 1]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.625 atol = (solver == :SCSMoment ? 1e-5 : 1e-7)
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2]
  PSD block sizes:
    [3 => 1]
> Clique #2: x[2], x[3]
  PSD block sizes:
    [3 => 1]"
    if optimize
        for solver in solvers
            @test poly_optimize(solver, sp).objective ≈ 0.625 atol = (solver == :SCSMoment ? 1e-4 : 1e-7)
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Example 3.4" begin
    DynamicPolynomials.@polyvar x[1:6]
    sp = Relaxation.SparsityCorrelativeTerm(
        poly_problem(1 + sum(x[i]^4 for i in 1:6) + x[1] * x[2] * x[3] + x[3] * x[4] * x[5] + x[3] * x[4] * x[6] +
                     x[3] * x[5] * x[6] + x[4] * x[5] * x[6]), 2
    )
    @test strRep(sp) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[3], x[4], x[5], x[6]
  PSD block sizes:
    [10 => 1, 5 => 1]
> Clique #2: x[1], x[2], x[3]
  PSD block sizes:
    [4 => 1, 2 => 3]"
    if optimize
        for solver in solvers
            @test poly_optimize(solver, sp).objective ≈ 0.5042475 atol = (solver == :SCSMoment ? 1e-3 : 1e-6)
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[3], x[4], x[5], x[6]
  PSD block sizes:
    [15 => 1]
> Clique #2: x[1], x[2], x[3]
  PSD block sizes:
    [4 => 1, 2 => 3]"
    if optimize
        for solver in solvers
            @test poly_optimize(solver, sp).objective ≈ 0.5042475 atol = (solver == :SCSMoment ? 1e-5 : 1e-7)
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[3], x[4], x[5], x[6]
  PSD block sizes:
    [15 => 1]
> Clique #2: x[1], x[2], x[3]
  PSD block sizes:
    [6 => 1, 2 => 2]"
    if optimize
        for solver in solvers
            @test poly_optimize(solver, sp).objective ≈ 0.5042475 atol = (solver == :SCSMoment ? 1e-6 : 1e-7)
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[3], x[4], x[5], x[6]
  PSD block sizes:
    [15 => 1]
> Clique #2: x[1], x[2], x[3]
  PSD block sizes:
    [6 => 1, 4 => 1]"
    if optimize
        for solver in solvers
            @test poly_optimize(solver, sp).objective ≈ 0.5042475 atol = (solver == :SCSMoment ? 1e-3 : 1e-6)
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Broyden Banded Function" begin
    n = 20
    DynamicPolynomials.@polyvar x[1:n]
    prob = poly_problem(sum((x[i] * (2 + 5x[i]^2) + 1 -
                             sum(j == i ? 0 : (1 + x[j]) * x[j] for j = max(1, i - 5):min(n, i + 1)))^2 for i = 1:n))
    sps = [Relaxation.SparsityCorrelativeTerm(prob, 3, chordal_completion=false, method=Relaxation.TERM_MODE_CLIQUES),
        Relaxation.SparsityCorrelativeTerm(prob, 3, chordal_completion=false, method=Relaxation.TERM_MODE_CHORDAL_CLIQUES),
        Relaxation.SparsityCorrelativeTerm(prob, 3)]
    @test strRep(sps[1]) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7]
  PSD block sizes:
    [10 => 14, 9 => 9, 8 => 10, 7 => 5, 6 => 2, 5 => 8, 4 => 14, 3 => 3, 1 => 15]
> Clique #2: x[2], x[3], x[4], x[5], x[6], x[7], x[8]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 9, 7 => 3, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #3: x[3], x[4], x[5], x[6], x[7], x[8], x[9]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #4: x[4], x[5], x[6], x[7], x[8], x[9], x[10]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #5: x[5], x[6], x[7], x[8], x[9], x[10], x[11]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #6: x[6], x[7], x[8], x[9], x[10], x[11], x[12]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #7: x[7], x[8], x[9], x[10], x[11], x[12], x[13]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #8: x[8], x[9], x[10], x[11], x[12], x[13], x[14]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #9: x[9], x[10], x[11], x[12], x[13], x[14], x[15]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #10: x[10], x[11], x[12], x[13], x[14], x[15], x[16]
  PSD block sizes:
    [10 => 14, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #11: x[11], x[12], x[13], x[14], x[15], x[16], x[17]
  PSD block sizes:
    [10 => 16, 9 => 7, 8 => 7, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #12: x[12], x[13], x[14], x[15], x[16], x[17], x[18]
  PSD block sizes:
    [10 => 8, 9 => 11, 8 => 4, 7 => 2, 5 => 6, 4 => 11, 3 => 3, 1 => 15]
> Clique #13: x[13], x[14], x[15], x[16], x[17], x[18], x[19]
  PSD block sizes:
    [10 => 6, 9 => 11, 8 => 4, 7 => 2, 6 => 6, 5 => 7, 4 => 10, 3 => 5, 1 => 15]
> Clique #14: x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [11 => 4, 10 => 5, 9 => 21, 8 => 9, 7 => 4, 6 => 6, 5 => 34, 4 => 34, 3 => 11, 1 => 35]"
    @test strRep(sps[2]) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7]
  PSD block sizes:
    [18 => 3, 17 => 3, 15 => 1, 14 => 2, 13 => 2, 12 => 2, 11 => 2, 10 => 3, 9 => 5, 7 => 1, 6 => 3, 4 => 2, 1 => 15]
> Clique #2: x[2], x[3], x[4], x[5], x[6], x[7], x[8]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 2, 13 => 3, 12 => 3, 11 => 3, 10 => 3, 9 => 6, 7 => 1, 6 => 3, 4 => 2, 1 => 15]
> Clique #3: x[3], x[4], x[5], x[6], x[7], x[8], x[9]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 2, 13 => 3, 12 => 3, 11 => 3, 10 => 3, 9 => 6, 7 => 1, 6 => 3, 4 => 2, 1 => 15]
> Clique #4: x[4], x[5], x[6], x[7], x[8], x[9], x[10]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 6, 7 => 1, 6 => 3, 4 => 2, 1 => 15]
> Clique #5: x[5], x[6], x[7], x[8], x[9], x[10], x[11]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 6, 7 => 1, 6 => 3, 4 => 2, 1 => 15]
> Clique #6: x[6], x[7], x[8], x[9], x[10], x[11], x[12]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 6, 8 => 1, 7 => 2, 6 => 3, 4 => 2, 1 => 15]
> Clique #7: x[7], x[8], x[9], x[10], x[11], x[12], x[13]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 6, 8 => 1, 7 => 2, 6 => 3, 4 => 2, 1 => 15]
> Clique #8: x[8], x[9], x[10], x[11], x[12], x[13], x[14]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 6, 8 => 1, 7 => 2, 6 => 3, 4 => 2, 1 => 15]
> Clique #9: x[9], x[10], x[11], x[12], x[13], x[14], x[15]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 6, 8 => 1, 7 => 2, 6 => 3, 4 => 2, 1 => 15]
> Clique #10: x[10], x[11], x[12], x[13], x[14], x[15], x[16]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 6, 8 => 1, 7 => 2, 6 => 3, 4 => 2, 1 => 15]
> Clique #11: x[11], x[12], x[13], x[14], x[15], x[16], x[17]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 6, 8 => 1, 7 => 2, 6 => 3, 4 => 2, 1 => 15]
> Clique #12: x[12], x[13], x[14], x[15], x[16], x[17], x[18]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 7, 8 => 1, 7 => 2, 6 => 3, 4 => 2, 1 => 15]
> Clique #13: x[13], x[14], x[15], x[16], x[17], x[18], x[19]
  PSD block sizes:
    [18 => 3, 17 => 4, 15 => 1, 14 => 3, 13 => 3, 12 => 4, 11 => 3, 10 => 3, 9 => 5, 8 => 1, 7 => 3, 6 => 3, 4 => 2, 1 => 15]
> Clique #14: x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [18 => 3, 17 => 5, 15 => 2, 14 => 6, 13 => 4, 12 => 7, 11 => 8, 10 => 5, 9 => 9, 8 => 7, 7 => 8, 6 => 11, 4 => 5, 1 => 35]"
    @test strRep(sps[3]) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #2: x[2], x[3], x[4], x[5], x[6], x[7], x[8]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #3: x[3], x[4], x[5], x[6], x[7], x[8], x[9]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #4: x[4], x[5], x[6], x[7], x[8], x[9], x[10]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #5: x[5], x[6], x[7], x[8], x[9], x[10], x[11]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #6: x[6], x[7], x[8], x[9], x[10], x[11], x[12]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #7: x[7], x[8], x[9], x[10], x[11], x[12], x[13]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #8: x[8], x[9], x[10], x[11], x[12], x[13], x[14]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #9: x[9], x[10], x[11], x[12], x[13], x[14], x[15]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #10: x[10], x[11], x[12], x[13], x[14], x[15], x[16]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #11: x[11], x[12], x[13], x[14], x[15], x[16], x[17]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #12: x[12], x[13], x[14], x[15], x[16], x[17], x[18]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #13: x[13], x[14], x[15], x[16], x[17], x[18], x[19]
  PSD block sizes:
    [85 => 1, 1 => 15]
> Clique #14: x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [85 => 1, 1 => 35]"

    # we test some of the optimziations, but some just take too long
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sps[1]).objective ≈ 0 atol = (solver ∈ (:HypatiaMoment, :SCSMoment) ? 3e-4 : 6e-6)
            end
        end

        for solver in solvers
            @testset let solver=solver
                @test(poly_optimize(solver, sps[2]).objective ≈ 0, atol = (occursin("Mosek", string(solver)) ? 2e-7 : 1e-4),
                    skip = (solver ∈ (:ClarabelMoment, :SCSMoment))) # Clarabel: 94 seconds; SCS: 85 seconds
            end
        end

        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, sps[3]).objective ≈ 0 atol = 2e-7
        :COPTMoment ∈ solvers && @test poly_optimize(:COPTMoment, sps[3]).objective ≈ 0 atol = 2e-7
        # MosekMoment: 77 seconds; HypatiaMoment, ClarabelMoment, SCSMoment: very slow (tested for Clarabel: 1h)
    end

    @test isnothing(iterate!(sps[1]))
    @test strRep(iterate!(sps[2])) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7]
  PSD block sizes:
    [36 => 1, 35 => 3, 34 => 1, 33 => 3, 32 => 3, 28 => 1, 27 => 1, 26 => 2, 25 => 1, 24 => 2, 23 => 1, 22 => 1, 21 => 1, 20 => 2, 19 => 4, 18 => 1, 17 => 1, 16 => 1, 14 => 2, 13 => 1, 12 => 2, 11 => 1, 10 => 1, 7 => 1, 6 => 3, 5 => 2, 4 => 1, 3 => 5, 1 => 4]
> Clique #2: x[2], x[3], x[4], x[5], x[6], x[7], x[8]
  PSD block sizes:
    [39 => 1, 38 => 1, 36 => 2, 35 => 1, 33 => 1, 30 => 1, 29 => 2, 27 => 1, 26 => 2, 25 => 2, 24 => 1, 23 => 1, 22 => 3, 21 => 3, 20 => 2, 19 => 1, 18 => 2, 17 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 1, 8 => 1, 6 => 4, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #3: x[3], x[4], x[5], x[6], x[7], x[8], x[9]
  PSD block sizes:
    [40 => 1, 38 => 2, 37 => 3, 35 => 4, 32 => 1, 31 => 1, 29 => 1, 28 => 1, 26 => 1, 25 => 1, 24 => 1, 23 => 5, 22 => 2, 20 => 2, 19 => 1, 18 => 1, 17 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 1, 8 => 1, 6 => 4, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #4: x[4], x[5], x[6], x[7], x[8], x[9], x[10]
  PSD block sizes:
    [42 => 1, 40 => 2, 39 => 3, 37 => 1, 35 => 2, 33 => 2, 31 => 1, 28 => 1, 26 => 3, 24 => 1, 23 => 4, 22 => 2, 20 => 1, 19 => 2, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 2, 8 => 1, 6 => 3, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #5: x[5], x[6], x[7], x[8], x[9], x[10], x[11]
  PSD block sizes:
    [41 => 1, 40 => 2, 39 => 1, 38 => 2, 37 => 1, 36 => 2, 35 => 2, 32 => 1, 31 => 1, 28 => 1, 27 => 1, 26 => 1, 24 => 1, 23 => 3, 22 => 1, 20 => 2, 19 => 1, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 2, 8 => 1, 6 => 3, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #6: x[6], x[7], x[8], x[9], x[10], x[11], x[12]
  PSD block sizes:
    [41 => 1, 40 => 2, 39 => 1, 38 => 2, 37 => 1, 36 => 2, 35 => 2, 32 => 1, 31 => 1, 28 => 1, 27 => 2, 26 => 1, 24 => 2, 23 => 3, 22 => 1, 20 => 2, 19 => 1, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 2, 8 => 1, 6 => 3, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #7: x[7], x[8], x[9], x[10], x[11], x[12], x[13]
  PSD block sizes:
    [41 => 1, 40 => 2, 39 => 1, 38 => 2, 37 => 1, 36 => 2, 35 => 2, 32 => 1, 31 => 1, 28 => 1, 27 => 2, 26 => 1, 24 => 2, 23 => 3, 22 => 1, 20 => 2, 19 => 1, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 2, 8 => 1, 6 => 3, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #8: x[8], x[9], x[10], x[11], x[12], x[13], x[14]
  PSD block sizes:
    [41 => 1, 40 => 2, 39 => 1, 38 => 2, 37 => 1, 36 => 2, 35 => 2, 32 => 1, 31 => 1, 28 => 1, 27 => 2, 26 => 1, 24 => 2, 23 => 3, 22 => 1, 20 => 2, 19 => 1, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 2, 8 => 1, 6 => 3, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #9: x[9], x[10], x[11], x[12], x[13], x[14], x[15]
  PSD block sizes:
    [41 => 1, 40 => 2, 39 => 1, 38 => 2, 37 => 1, 36 => 2, 35 => 2, 32 => 1, 31 => 1, 28 => 1, 27 => 2, 26 => 1, 24 => 2, 23 => 3, 22 => 1, 20 => 2, 19 => 1, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 2, 8 => 1, 6 => 3, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #10: x[10], x[11], x[12], x[13], x[14], x[15], x[16]
  PSD block sizes:
    [41 => 1, 40 => 2, 39 => 1, 38 => 2, 37 => 1, 36 => 2, 35 => 2, 32 => 1, 31 => 1, 28 => 1, 27 => 2, 26 => 1, 24 => 2, 23 => 3, 22 => 1, 20 => 2, 19 => 1, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 2, 8 => 1, 6 => 3, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #11: x[11], x[12], x[13], x[14], x[15], x[16], x[17]
  PSD block sizes:
    [41 => 1, 40 => 2, 39 => 1, 38 => 2, 37 => 1, 36 => 2, 35 => 2, 32 => 1, 31 => 1, 28 => 1, 27 => 2, 26 => 1, 24 => 2, 23 => 3, 22 => 1, 20 => 2, 19 => 1, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 11 => 1, 10 => 2, 8 => 1, 6 => 3, 5 => 2, 4 => 3, 3 => 2, 1 => 4]
> Clique #12: x[12], x[13], x[14], x[15], x[16], x[17], x[18]
  PSD block sizes:
    [41 => 1, 40 => 2, 39 => 2, 37 => 1, 36 => 3, 35 => 1, 32 => 1, 31 => 1, 28 => 2, 27 => 2, 26 => 1, 25 => 1, 24 => 1, 23 => 3, 22 => 2, 20 => 2, 19 => 1, 18 => 1, 16 => 1, 15 => 1, 14 => 1, 12 => 1, 11 => 1, 10 => 1, 8 => 1, 6 => 4, 5 => 1, 4 => 3, 3 => 2, 1 => 4]
> Clique #13: x[13], x[14], x[15], x[16], x[17], x[18], x[19]
  PSD block sizes:
    [40 => 2, 39 => 2, 38 => 1, 37 => 1, 36 => 1, 35 => 1, 34 => 1, 33 => 2, 31 => 1, 28 => 2, 27 => 1, 26 => 1, 25 => 1, 24 => 3, 23 => 2, 22 => 3, 20 => 2, 19 => 3, 18 => 2, 16 => 1, 15 => 1, 14 => 1, 12 => 1, 11 => 1, 10 => 1, 8 => 1, 6 => 4, 5 => 1, 4 => 3, 3 => 2, 1 => 4]
> Clique #14: x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [37 => 4, 36 => 1, 35 => 3, 34 => 1, 33 => 3, 31 => 1, 29 => 2, 28 => 1, 27 => 1, 26 => 1, 25 => 2, 24 => 4, 23 => 3, 21 => 2, 19 => 2, 18 => 3, 17 => 1, 16 => 2, 15 => 2, 14 => 2, 13 => 1, 11 => 1, 10 => 3, 9 => 1, 8 => 2, 7 => 4, 6 => 3, 5 => 6, 4 => 3, 3 => 11, 1 => 4]"
    @test strRep(iterate!(sps[3])) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7]
  PSD block sizes:
    [120 => 1]
> Clique #2: x[2], x[3], x[4], x[5], x[6], x[7], x[8]
  PSD block sizes:
    [120 => 1]
> Clique #3: x[3], x[4], x[5], x[6], x[7], x[8], x[9]
  PSD block sizes:
    [120 => 1]
> Clique #4: x[4], x[5], x[6], x[7], x[8], x[9], x[10]
  PSD block sizes:
    [120 => 1]
> Clique #5: x[5], x[6], x[7], x[8], x[9], x[10], x[11]
  PSD block sizes:
    [120 => 1]
> Clique #6: x[6], x[7], x[8], x[9], x[10], x[11], x[12]
  PSD block sizes:
    [120 => 1]
> Clique #7: x[7], x[8], x[9], x[10], x[11], x[12], x[13]
  PSD block sizes:
    [120 => 1]
> Clique #8: x[8], x[9], x[10], x[11], x[12], x[13], x[14]
  PSD block sizes:
    [120 => 1]
> Clique #9: x[9], x[10], x[11], x[12], x[13], x[14], x[15]
  PSD block sizes:
    [120 => 1]
> Clique #10: x[10], x[11], x[12], x[13], x[14], x[15], x[16]
  PSD block sizes:
    [120 => 1]
> Clique #11: x[11], x[12], x[13], x[14], x[15], x[16], x[17]
  PSD block sizes:
    [120 => 1]
> Clique #12: x[12], x[13], x[14], x[15], x[16], x[17], x[18]
  PSD block sizes:
    [120 => 1]
> Clique #13: x[13], x[14], x[15], x[16], x[17], x[18], x[19]
  PSD block sizes:
    [120 => 1]
> Clique #14: x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [120 => 1]"
end

@testset "Generalized Rosenbrock" begin
    n = 40
    DynamicPolynomials.@polyvar x[1:n]
    sp = Relaxation.SparsityCorrelativeTerm(
        poly_problem(1 + sum(100 * (x[i] - x[i-1]^2)^2 + (1 - x[i])^2 for i in 2:n),
                     nonneg=[1 - sum(x[i]^2 for i in (20j-19):(20j)) for j in 1:Int(n / 20)]),
        method=Relaxation.TERM_MODE_CHORDAL_CLIQUES
    )
    @test strRep(sp) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [21 => 1, 3 => 19, 2 => 37, 1 => 172]
> Clique #2: x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31], x[32], x[33], x[34], x[35], x[36], x[37], x[38], x[39], x[40]
  PSD block sizes:
    [21 => 1, 3 => 19, 2 => 38, 1 => 171]
> Clique #3: x[20], x[21]
  PSD block sizes:
    [3 => 2, 2 => 3]"

    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 38.0494 atol = solver==:SCSMoment ? 2e-1 : 1e-4
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [40 => 1, 3 => 171, 2 => 37, 1 => 1]
> Clique #2: x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31], x[32], x[33], x[34], x[35], x[36], x[37], x[38], x[39], x[40]
  PSD block sizes:
    [41 => 1, 3 => 190, 2 => 19]
> Clique #3: x[20], x[21]
  PSD block sizes:
    [4 => 1, 2 => 3]"

    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 38.0514 atol = solver==:SCSMoment ? 1e-2 : 2e-4
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [40 => 1, 20 => 1, 4 => 171, 2 => 19, 1 => 1]
> Clique #2: x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31], x[32], x[33], x[34], x[35], x[36], x[37], x[38], x[39], x[40]
  PSD block sizes:
    [41 => 1, 21 => 1, 4 => 190]
> Clique #3: x[20], x[21]
  PSD block sizes:
    [4 => 1, 2 => 1]"

    if optimize
        for solver in solvers
            @testset let solver=solver
                @test(poly_optimize(solver, sp).objective ≈ 38.0514,
                    atol = (solver==:SCSMoment ? 2e-2 : 2e-4), skip = solver==:COPTMoment)
            end
        end
        # COPT: 170s
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [149 => 1, 133 => 2, 89 => 5, 58 => 10, 20 => 2, 1 => 1]
> Clique #2: x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31], x[32], x[33], x[34], x[35], x[36], x[37], x[38], x[39], x[40]
  PSD block sizes:
    [165 => 1, 141 => 2, 107 => 1, 93 => 4, 77 => 1, 60 => 10, 21 => 1]
> Clique #3: x[20], x[21]
  PSD block sizes:
    [4 => 1, 2 => 1]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20]
  PSD block sizes:
    [211 => 1, 20 => 2, 1 => 1]
> Clique #2: x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], x[29], x[30], x[31], x[32], x[33], x[34], x[35], x[36], x[37], x[38], x[39], x[40]
  PSD block sizes:
    [231 => 1, 21 => 1]
> Clique #3: x[20], x[21]
  PSD block sizes:
    [4 => 1, 2 => 1]"

    @test isnothing(iterate!(sp))
end

@testset "Example 4.4 from Zhen, Fantuzzi, Papachristodoulou review" begin
    DynamicPolynomials.@polyvar x[1:4]
    sp = Relaxation.SparsityCorrelativeTerm(
        poly_problem(2 + x[1]^2 * x[4]^2 * (x[1]^2 * x[4]^2 - 1) - x[1]^2 + x[1]^4 +
                     sum((x[i]^2 * x[i-1]^2 * (x[i]^2 * x[i-1]^2 - 1) - x[i]^2 + x[i]^4 for i in 2:4))), 4
    )
    @test strRep(sp) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2], x[4]
  PSD block sizes:
    [10 => 1, 4 => 6, 1 => 1]
> Clique #2: x[2], x[3], x[4]
  PSD block sizes:
    [10 => 1, 4 => 6, 1 => 1]"

    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.110118 atol = solver==:SCSMoment ? 1e-4 : 1e-6
            end
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Something with matrices" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityCorrelativeTerm(poly_problem(1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1] * x[2] + x[2],
                                                         psd=[[x[1] x[2]; x[2] x[1]]]), 2)
    @test strRep(sp) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: x[1], x[2]
  PSD block sizes:
    [6 => 1, 3 => 1]
> Clique #2: x[3]
  PSD block sizes:
    [2 => 1, 1 => 1]"

    if optimize
        for solver in solvers
            @testset let solver=solver
                @test(poly_optimize(solver, sp; (solver == :HypatiaMoment ? Dict(:dense => true) : Dict())...).objective ≈ 0.283871,
                    atol = solver==:SCSMoment ? 1e-4 : 1e-6)
            end
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Complex" begin
    DynamicPolynomials.@complex_polyvar z[1:2]
    sp = Relaxation.SparsityCorrelativeTerm(poly_problem(z[1] + conj(z[1]),
                                                         nonneg=[1 - z[1] * conj(z[1]) - z[2] * conj(z[2])]), 2)
    @test strRep(sp) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: z[1], z[2]
  PSD block sizes:
    [2 => 2, 1 => 5]"

    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ -2 atol = solver==:SCSMoment ? 1e-4 : 1e-6
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityCorrelativeTerm of a polynomial optimization problem
> Clique #1: z[1], z[2]
  PSD block sizes:
    [3 => 1, 2 => 2, 1 => 2]"

    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ -2 atol = solver==:SCSMoment ? 1e-4 : 1e-7
            end
        end
    end

    @test isnothing(iterate!(sp))
end
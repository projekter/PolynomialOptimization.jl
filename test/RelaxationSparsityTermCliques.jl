include("./shared.jl")

# These tests are a bit problematic, since the chordal extension is not unique. Even if the algorithm is fully specified, the
# actual results depend on the internal structure of the graphs (i.e., term order). So do not expect identical results when
# comparing with differently ordered terms.

@testset "Example 3.6" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityTermChordal(Relaxation.Custom(
        poly_problem(x[1]^2 - 2x[1] * x[2] + 3x[2]^2 - 2x[1]^2 * x[2] + 2x[1]^2 * x[2]^2 - 2x[2] * x[3] +
                     6x[3]^2 + 18x[2]^2 * x[3] - 54x[2] * x[3]^2 + 142x[2]^2 * x[3]^2),
        monomial_type(x[1])[1, x[1], x[2], x[3], x[1] * x[2], x[2] * x[3]])
    )
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [3 => 4]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ -0.0035512 atol = solver==:SCSMoment ? 1e-4 : 1e-6
            end
        end
    end

    # paper says that the iterations terminate. However, this is a result of the paper always considering the constant term to
    # be a part of the objective. We don't do this, so we can go to the next level.
    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [4 => 2, 3 => 1]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0 atol = solver==:SCSMoment ? 1e-6 : 5e-8
            end
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Example 7.2 from TSSOS with block completion" begin
    DynamicPolynomials.@polyvar x[1:3] y[1:3]
    sp = Relaxation.SparsityTermChordal(
        poly_problem(27 - ((x[1] - x[2])^2 + (y[1] - y[2])^2) * ((x[1] - x[3])^2 + (y[1] - y[3])^2) *
                          ((x[2] - x[3])^2 + (y[2] - y[3])^2),
                     zero=[x[1]^2 + y[1]^2 + x[2]^2 + y[2]^2 + x[3]^2 + y[3]^2 - 3]),
        3
    )
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6]
PSD block sizes:
  [24 => 1, 23 => 2, 21 => 1, 20 => 2, 18 => 2, 7 => 7, 1 => 15]
Free block sizes:
  [12 => 1, 9 => 1, 7 => 1, 1 => 6]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test(poly_optimize(solver, sp;
                    (solver == :HypatiaMoment ? Dict(:preprocess => true) : Dict())...).objective ≈ 0,
                    atol=solver==:SCSMoment ? 5e-5 : 9e-7)
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6]
PSD block sizes:
  [29 => 6, 12 => 1, 9 => 1, 7 => 1]
Free block sizes:
  [13 => 1, 9 => 1, 3 => 2]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6]
PSD block sizes:
  [31 => 2, 13 => 1, 9 => 1]
Free block sizes:
  [13 => 1, 9 => 1, 3 => 2]"

    @test isnothing(iterate!(sp))
end

@testset "Example 4.5 from Zhen, Fantuzzi, Papachristodoulou review" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityTermChordal(poly_problem(1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]^2 * x[2]^2 + x[1]^2 * x[3]^2 +
                                                    x[2]^2 * x[3]^2 + x[2] * x[3]), 2)
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [4 => 1, 2 => 2, 1 => 3]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.91666667 atol = 1e-6
            end
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Something with matrices" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityTermChordal(poly_problem(1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1] * x[2] * x[3] + x[2],
                                                    psd=[[x[1] x[2]; x[2] x[1]]]), 2)
    @test strRep(groupings(sp)) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[x₁, x₂, x₃]

Block groupings
===============
Objective: 6 blocks
  4 [1, x₃², x₂², x₁²]
  2 [1, x₂]
  2 [1, x₁]
  2 [x₃, x₁x₂]
  2 [x₂, x₁x₃]
  2 [x₁, x₂x₃]
Semidefinite constraint #1: 2 blocks
  3 [1, x₂, x₁]
  3 [x₃, x₂, x₁]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.5355788 atol = solver==:SCSMoment ? 1e-4 : 1e-7
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [5 => 2, 4 => 2, 3 => 4]"

    @test isnothing(iterate!(sp))
end

@testset "Other example" begin
    n = 7
    DynamicPolynomials.@polyvar x[1:n]
    sp = Relaxation.SparsityTermChordal(
        poly_problem(
            sum((x[i] * (2 + 5x[i]^2) + 1 -
                sum((j != i) && (max(1, i - 5) ≤ j ≤ min(n, i + 1)) ? (1 + x[j]) * x[j] : 0
                    for j = 1:n))^2 for i = 1:n)
        ),
        3
    )
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [18 => 3, 17 => 4, 15 => 2, 14 => 3, 13 => 3, 12 => 4, 11 => 5, 10 => 4, 9 => 8, 8 => 3, 7 => 4, 6 => 11, 4 => 5, 1 => 35]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0 atol = 2e-5
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [33 => 3, 32 => 2, 31 => 2, 30 => 1, 28 => 1, 27 => 1, 26 => 1, 25 => 2, 24 => 2, 23 => 5, 22 => 2, 21 => 3, 20 => 1, 19 => 3, 18 => 4, 17 => 1, 16 => 4, 15 => 1, 14 => 3, 13 => 2, 12 => 2, 11 => 1, 10 => 1, 7 => 5, 6 => 4, 5 => 8, 4 => 1, 3 => 14, 1 => 4]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [43 => 1, 41 => 3, 40 => 1, 38 => 1, 36 => 1, 35 => 2, 34 => 2, 33 => 2, 32 => 4, 30 => 1, 29 => 4, 28 => 2, 27 => 3, 26 => 3, 25 => 1, 24 => 3, 23 => 2, 22 => 4, 21 => 4, 20 => 1, 19 => 2, 18 => 4, 17 => 5, 16 => 4, 15 => 5, 14 => 4, 12 => 1, 10 => 1, 9 => 3]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [49 => 3, 48 => 2, 46 => 1, 44 => 3, 43 => 1, 39 => 2, 38 => 2, 36 => 3, 34 => 5, 33 => 3, 32 => 2, 31 => 3, 30 => 1, 29 => 4, 28 => 4, 27 => 3, 26 => 1, 25 => 4, 24 => 4, 23 => 5, 22 => 3, 21 => 4, 20 => 3]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [50 => 1, 49 => 3, 48 => 1, 45 => 4, 40 => 2, 39 => 1, 38 => 2, 37 => 6, 36 => 2, 35 => 6, 34 => 3, 33 => 4, 32 => 6, 31 => 3, 30 => 8, 29 => 1, 27 => 1, 26 => 1, 25 => 4, 24 => 2, 23 => 1, 22 => 1]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [50 => 1, 49 => 3, 48 => 1, 45 => 4, 44 => 3, 43 => 2, 40 => 6, 39 => 4, 38 => 3, 37 => 3, 36 => 2, 35 => 6, 34 => 3, 33 => 5, 32 => 1, 31 => 1, 30 => 5, 27 => 4, 26 => 4, 25 => 1, 24 => 1]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [51 => 1, 50 => 1, 49 => 1, 48 => 1, 46 => 1, 45 => 5, 44 => 1, 43 => 1, 42 => 1, 41 => 3, 40 => 5, 39 => 4, 37 => 4, 36 => 1, 35 => 8, 34 => 2, 33 => 3, 32 => 1, 31 => 2, 30 => 3, 29 => 1, 28 => 3, 27 => 2, 26 => 2, 25 => 2]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [51 => 1, 50 => 1, 49 => 1, 48 => 1, 46 => 1, 45 => 5, 44 => 2, 42 => 1, 41 => 4, 40 => 5, 39 => 3, 37 => 4, 36 => 1, 35 => 8, 34 => 3, 33 => 2, 32 => 1, 31 => 3, 30 => 2, 29 => 1, 28 => 3, 27 => 3, 26 => 1, 25 => 2]"

    @test isnothing(iterate!(sp))
end

@testset "Complex-valued" begin
    DynamicPolynomials.@complex_polyvar z[1:2]
    sp = Relaxation.SparsityTermChordal(poly_problem(z[1] + conj(z[1]), nonneg=[1 - z[1] * conj(z[1]) - z[2] * conj(z[2])]), 2)
    @test strRep(groupings(sp)) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[z₁, z₂]

Block groupings
===============
Objective: 5 blocks
  2 [1, z₁]
  1 [z₂]
  1 [z₂²]
  1 [z₁z₂]
  1 [z₁²]
Nonnegative constraint #1: 2 blocks
  2 [1, z₁]
  1 [z₂]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ -2 atol = solver==:SCSMoment ? 1e-4 : 1e-6
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  z[1], z[2]
PSD block sizes:
  [2 => 4, 1 => 2]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ -2 atol = 1e-4
            end
        end
    end

    @test isnothing(iterate!(sp))
end
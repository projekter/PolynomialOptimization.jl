include("./shared.jl")

filter!(s -> !occursin("SCS", string(s)), solvers) # precision issues

@testset "Example 4.3" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityTermBlock(poly_problem(1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1] * x[2] * x[3] + x[2]), 2)
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [6 => 1, 2 => 2]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.4752747 atol = 1e-6
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [6 => 1, 4 => 1]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.4752747 atol = 1e-6
            end
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Example 4.5 from Zhen, Fantuzzi, Papachristodoulou review" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityTermBlock(poly_problem(1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]^2 * x[2]^2 + x[1]^2 * x[3]^2 +
                                                  x[2]^2 * x[3]^2 + x[2] * x[3]), 2)
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [5 => 1, 2 => 1, 1 => 3]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.91666667 atol = 1e-6
            end
        end
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [5 => 1, 2 => 2, 1 => 1]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.91666667 atol = 1e-6
            end
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Example 7.2" begin
    DynamicPolynomials.@polyvar x[1:3] y[1:3]
    sp = Relaxation.SparsityTermBlock(
      poly_problem(27 - ((x[1] - x[2])^2 + (y[1] - y[2])^2) * ((x[1] - x[3])^2 + (y[1] - y[3])^2) *
                        ((x[2] - x[3])^2 + (y[2] - y[3])^2),
                   zero=[x[1]^2 + y[1]^2 + x[2]^2 + y[2]^2 + x[3]^2 + y[3]^2 - 3]),
      3
    )
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6]
PSD block sizes:
  [31 => 2, 7 => 1, 1 => 15]
Free block sizes:
  [13 => 1, 9 => 1, 1 => 6]"
    if optimize
        :MosekMoment ∈ solvers && @test poly_optimize(:MosekMoment, sp).objective ≈ 0 atol = 1e-9
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, sp).objective ≈ 0 atol = 1e-7
        :COSMOMoment ∈ solvers && @test poly_optimize(:COSMOMoment, sp).objective ≈ 0 atol = 1e-6
        :HypatiaMoment ∈ solvers && @test poly_optimize(:HypatiaMoment, sp).objective ≈ 0 atol = 2e-6
        :COPTSOS ∈ solvers && @test poly_optimize(:COPTSOS, sp).objective ≈ 0 atol = 1e-6
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6]
PSD block sizes:
  [31 => 2, 13 => 1, 9 => 1]
Free block sizes:
  [13 => 1, 9 => 1, 3 => 2]"

    @test isnothing(iterate!(sp))
end

@testset "Something with matrices" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityTermBlock(poly_problem(1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1] * x[2] * x[3] + x[2],
                                                  psd=[[x[1] x[2]; x[2] x[1]]]), 2)
    @test strRep(groupings(sp)) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[x₁, x₂, x₃]

Block groupings
===============
Objective: 2 blocks
  8 [1, x₂, x₁, x₃², x₂x₃, x₂², x₁x₃, x₁²]
  2 [x₃, x₁x₂]
Semidefinite constraint #1: 1 block
  4 [1, x₃, x₂, x₁]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.5355788 atol = 1e-6
            end
        end
    end

    @test strRep(groupings(iterate!(sp))) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[x₁, x₂, x₃]

Block groupings
===============
Objective: 1 block
  10 [1, x₃, x₂, x₁, x₃², x₂x₃, x₂², x₁x₃, x₁x₂, x₁²]
Semidefinite constraint #1: 1 block
  4 [1, x₃, x₂, x₁]"
    @test isnothing(iterate!(sp))
end

@testset "Other example" begin
    n = 7
    DynamicPolynomials.@polyvar x[1:n]
    sp = Relaxation.SparsityTermBlock(poly_problem(
        sum((x[i] * (2 + 5x[i]^2) + 1 -
             sum((j != i) && (max(1, i - 5) ≤ j ≤ min(n, i + 1)) ? (1 + x[j]) * x[j] : 0
                 for j = 1:n))^2 for i = 1:n),
    ), 3)
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [85 => 1, 1 => 35]"
    if optimize
        :MosekMoment ∈ solvers && @test poly_optimize(:MosekMoment, sp).objective ≈ 0 atol = 2e-8
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, sp).objective ≈ 0 atol = 1e-8
        :COSMOMoment ∈ solvers && @test poly_optimize(:COSMOMoment, sp).objective ≈ 0 atol = 1e-4
        :HypatiaMoment ∈ solvers && @test poly_optimize(:HypatiaMoment, sp, dense=true).objective ≈ 0 atol = 1e-7
        :COPTSOS ∈ solvers && @test poly_optimize(:COPTSOS, sp).objective ≈ 0 atol = 2e-7
    end

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [120 => 1]"
end

@testset "Complex-valued" begin
    DynamicPolynomials.@complex_polyvar z[1:2]
    sp = Relaxation.SparsityTermBlock(poly_problem(z[1] + conj(z[1]), nonneg=[1 - z[1] * conj(z[1]) - z[2] * conj(z[2])]), 2)
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
        :MosekMoment ∈ solvers && @test poly_optimize(:MosekMoment, sp).objective ≈ -2 atol = 1e-8
        :HypatiaMoment ∈ solvers && @test poly_optimize(:HypatiaMoment, sp, dense=true).objective ≈ -2 atol = 1e-8
    end

    @test strRep(groupings(iterate!(sp))) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[z₁, z₂]

Block groupings
===============
Objective: 3 blocks
  3 [1, z₁, z₁²]
  2 [z₂, z₁z₂]
  1 [z₂²]
Nonnegative constraint #1: 2 blocks
  2 [1, z₁]
  1 [z₂]"
    if optimize
        :MosekMoment ∈ solvers && @test poly_optimize(:MosekMoment, sp).objective ≈ -2 atol = 1e-7
        :HypatiaMoment ∈ solvers && @test poly_optimize(:HypatiaMoment, sp, dense=true).objective ≈ -2 atol = 1e-7
    end

    @test isnothing(iterate!(sp))
end
include("./shared.jl")

@testset "Example 4.3" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityTermBlock(poly_problem(1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1] * x[2] * x[3] + x[2]), 2)
    @test strRep(sp) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [6 => 1, 2 => 2]"
    @test poly_optimize(:Clarabel, sp).objective ≈ 0.4752747 atol = 1e-6

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [6 => 1, 4 => 1]"
    @test poly_optimize(:Clarabel, sp).objective ≈ 0.4752747 atol = 1e-6

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
    @test poly_optimize(:Clarabel, sp).objective ≈ 0.91666667 atol = 1e-6

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [5 => 1, 2 => 2, 1 => 1]"
    @test poly_optimize(:Clarabel, sp).objective ≈ 0.91666667 atol = 1e-6

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
    @test poly_optimize(:COPT, sp).objective ≈ 0 atol = 1e-6

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
  - 2 indices: 4 [1, x₃, x₂, x₁]"
    @test poly_optimize(:COPT, sp).objective ≈ 0.5355788 atol = 1e-7

    @test strRep(groupings(iterate!(sp))) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[x₁, x₂, x₃]

Block groupings
===============
Objective: 1 block
  10 [1, x₃, x₂, x₁, x₃², x₂x₃, x₂², x₁x₃, x₁x₂, x₁²]
Semidefinite constraint #1: 1 block
  - 2 indices: 4 [1, x₃, x₂, x₁]"
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
    @test poly_optimize(:MosekMoment, sp).objective ≈ 0 atol = 2e-8 skip = !have_mosek
    @test poly_optimize(:MosekSOS, sp).objective ≈ 0 atol = 1e-8 skip = !have_mosek
    @test poly_optimize(:LoRADS, sp).objective ≈ 0 atol = 5e-5 skip = !have_lorads
    @test poly_optimize(:COPT, sp).objective ≈ 0 atol = 2e-7

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
    @test poly_optimize(:Clarabel, sp).objective ≈ -2 atol = 2e-8
    @test poly_optimize(:Hypatia, sp, dense=true).objective ≈ -2 atol = 1e-8

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
    @test poly_optimize(:Clarabel, sp).objective ≈ -2 atol = 1e-8
    @test poly_optimize(:HypatiaMoment, sp, dense=true).objective ≈ -2 atol = 1e-7

    @test isnothing(iterate!(sp))
end

@testset "Matrix term sparsity" begin
    # There's no usable example in the term sparsity paper, as they all use a matrix-valued objective, which is not supported
    # by PolynomialOptimization. We take a simplified version of their first example in 7.1.
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityTermBlock(poly_problem(zero(polynomial_type(x[1])),
                                                   psd=[[1-x[1]^2-x[2]^2 x[2]*x[3]; x[2]*x[3] 1-x[3]^2]]), 2,
                                      constraints_in_support=false) # for compatibility with the paper and TSSOS
    @test strRep(groupings(sp)) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[x₁, x₂, x₃]

Block groupings
===============
Objective: 7 blocks
  4 [1, x₃², x₂², x₁²]
  1 [x₃]
  1 [x₂]
  1 [x₁]
  1 [x₂x₃]
  1 [x₁x₃]
  1 [x₁x₂]
Semidefinite constraint #1: 6 blocks
  - index 1: 1 [x₃]
    index 2: 1 [x₂]
  - index 1: 1 [x₂]
    index 2: 1 [x₃]
  - index 2: 1 [1]
  - index 2: 1 [x₁]
  - index 1: 1 [1]
  - index 1: 1 [x₁]"
    # Note: while this is cross-verified with TSSOS, I don't trust the implementation there - it discards diagonal entries when
    # constructing the graph. In this case, no diagonal entry would make it anyway, but this will not always be the case.
    @test poly_optimize(:Clarabel, sp).objective ≈ 0 atol=1e-10
    @test isnothing(iterate!(sp))

    sp = Relaxation.SparsityTermBlock(poly_problem(sp), 2) # now with constraints considered in the support
    @test strRep(groupings(sp)) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[x₁, x₂, x₃]

Block groupings
===============
Objective: 5 blocks
  5 [1, x₃², x₂x₃, x₂², x₁²]
  2 [x₃, x₂]
  1 [x₁]
  1 [x₁x₃]
  1 [x₁x₂]
Semidefinite constraint #1: 4 blocks
  - 2 indices: 2 [x₃, x₂]
  - 2 indices: 1 [1]
  - index 2: 1 [x₁]
  - index 1: 1 [x₁]"
    @test poly_optimize(:Clarabel, sp).objective ≈ 0 atol=1e-10

    @test strRep(groupings(iterate!(sp))) == "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[x₁, x₂, x₃]

Block groupings
===============
Objective: 4 blocks
  5 [1, x₃², x₂x₃, x₂², x₁²]
  2 [x₃, x₂]
  2 [x₁x₃, x₁x₂]
  1 [x₁]
Semidefinite constraint #1: 3 blocks
  - 2 indices: 2 [x₃, x₂]
  - 2 indices: 1 [1]
  - 2 indices: 1 [x₁]"
    @test poly_optimize(:Clarabel, sp).objective ≈ 0 atol=1e-10

    @test isnothing(iterate!(sp))
end
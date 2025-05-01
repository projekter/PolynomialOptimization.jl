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
    @test poly_optimize(:Clarabel, sp).objective ≈ -0.0035512 atol = 1e-6

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
  [23 => 4, 20 => 2, 18 => 2, 7 => 7, 1 => 15]
Free block sizes:
  [12 => 1, 9 => 1, 7 => 1, 1 => 6]"
    @test poly_optimize(:Clarabel, sp).objective ≈ 0 atol = 5e-6

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
    @test poly_optimize(:Clarabel, sp).objective ≈ 0.916666667 atol = 1e-9

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
    @test poly_optimize(:Clarabel, sp).objective ≈ 0.5355788 atol = 1e-7

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
PSD block sizes:
  [6 => 2, 5 => 2, 4 => 2, 3 => 2]"

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
  [18 => 3, 17 => 5, 15 => 1, 14 => 2, 13 => 4, 12 => 5, 11 => 4, 10 => 4, 9 => 8, 8 => 3, 7 => 4, 6 => 11, 4 => 5, 1 => 35]"
    @test poly_optimize(:Clarabel, sp).objective ≈ 0 atol = 1e-5

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [33 => 2, 32 => 2, 31 => 1, 30 => 2, 28 => 1, 26 => 3, 25 => 3, 24 => 1, 23 => 1, 22 => 3, 21 => 2, 20 => 5, 19 => 2, 18 => 1, 17 => 2, 16 => 4, 15 => 3, 14 => 1, 13 => 2, 12 => 2, 11 => 1, 10 => 1, 7 => 5, 6 => 5, 5 => 8, 3 => 14, 1 => 4]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [44 => 1, 43 => 2, 41 => 1, 40 => 2, 39 => 1, 37 => 1, 36 => 1, 34 => 2, 33 => 1, 32 => 3, 31 => 4, 30 => 5, 29 => 2, 28 => 4, 26 => 3, 25 => 1, 24 => 5, 22 => 1, 21 => 1, 20 => 3, 19 => 1, 18 => 3, 17 => 5, 16 => 10, 15 => 2, 14 => 1, 13 => 2, 12 => 1, 11 => 2, 10 => 1, 9 => 1]"

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3], x[4], x[5], x[6], x[7]
PSD block sizes:
  [45 => 4, 44 => 3, 43 => 8, 42 => 2, 41 => 3, 39 => 3, 37 => 3, 36 => 2, 35 => 5, 33 => 1, 32 => 1, 31 => 7, 30 => 4, 29 => 4, 28 => 4, 27 => 3, 26 => 3, 25 => 2, 24 => 4, 22 => 2, 21 => 2, 18 => 1]"

  # this will go on for some time
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
    @test poly_optimize(:Clarabel, sp).objective ≈ -2 atol = 1e-7

    @test strRep(iterate!(sp)) == "Relaxation.SparsityTerm of a polynomial optimization problem
Variable cliques:
  z[1], z[2]
PSD block sizes:
  [2 => 4, 1 => 2]"
    @test poly_optimize(:Clarabel, sp).objective ≈ -2 atol = 1e-10

    @test isnothing(iterate!(sp))
end
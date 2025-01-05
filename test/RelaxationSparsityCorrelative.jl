include("./shared.jl")

# filter out very slow solvers. Indeed Mosek in its moment formulation sucks.
filter!(s -> s != :MosekMoment && s != :HypatiaMoment, solvers)

@testset "Example 6.1 chained singular from correlative sparsity paper" begin
    for n in (16, 40, 100, 200, 400)
        DynamicPolynomials.@polyvar x[1:n]
        sp = Relaxation.SparsityCorrelative(poly_problem(sum((x[i] + 10x[i+1])^2 + 5(x[i+2] - x[i+3])^2 +
                                                             (x[i+1] - 2x[i+2])^4 + 10(x[i] - 10x[i+3])^4 for i in 1:2:n-3)))
        @test StatsBase.countmap(length.(groupings(sp).var_cliques)) == Dict(3 => n -2)
        if optimize
            for solver in solvers
                # all moment-based solvers fail miserably on this problem. However, they will all report the issues.
                @testset let n=n, solver=solver
                    @test poly_optimize(solver, sp).objective ≈ 0 atol = 6e-5 skip=occursin("Moment", string(solver))
                end
            end
        end
    end
end

@testset "Example 6.1 Broyden banded function from correlative sparsity paper" begin
    for (n, cl) in ((6, Dict(6 => 1)), (7, Dict(7 => 1)), (8, Dict(7 => 2)), (9, Dict(7 => 3)), (10, Dict(7 => 4)))
        DynamicPolynomials.@polyvar x[1:n]
        sp = Relaxation.SparsityCorrelative(
            poly_problem(sum((x[i] * (2 + 5x[i]^2) + 1 -
                              sum(j == i ? 0 : (1 + x[j]) * x[j] for j in max(1, i - 5):min(n, i + 1)))^2 for i = 1:n))
        )
        @test StatsBase.countmap(length.(groupings(sp).var_cliques)) == cl
        if optimize
            for solver in solvers
                @testset let n=n, solver=solver
                    @test poly_optimize(solver, sp).objective ≈ 0 atol = 2e-5 skip=solver==:ClarabelMoment
                end
            end
        end
    end
end

@testset "Example 6.1 Broyden tridiagonal function from correlative sparsity paper" begin
    for n in 600:100:1000
        DynamicPolynomials.@polyvar x[1:n]
        sp = Relaxation.SparsityCorrelative(
            poly_problem(((3 - 2x[1]) * x[1] - 2x[2] + 1)^2 +
                sum(((3 - 2x[i]) * x[i] - x[i-1] - 2x[i+1] + 1)^2 for i in 2:n-1) + ((3 - 2x[n]) * x[n] - x[n-1] + 1)^2)
        )
        @test StatsBase.countmap(length.(groupings(sp).var_cliques)) == Dict(3 => n -2)
        if optimize
            for solver in solvers
                @testset let n=n, solver=solver
                    @test poly_optimize(solver, sp).objective ≈ 0 atol = 5e-5 skip=solver==:SCSMoment
                end
            end
        end
    end
end

@testset "Example 6.1 Chained Wood function from correlative sparsity paper" begin
    for n in 600:100:1000
        DynamicPolynomials.@polyvar x[1:n]
        sp = Relaxation.SparsityCorrelative(
            poly_problem(1 + sum(100(x[i+1] - x[i]^2)^2 + (1 - x[i])^2 + 90(x[i+3] - x[i+2]^2)^2 +
                (1 - x[i+2])^2 + 10(x[i+1] + x[i+3] - 2)^2 + 0.1(x[i+1] - x[i+3])^2 for i in 1:2:n-3))
        )
        @test StatsBase.countmap(length.(groupings(sp).var_cliques)) == Dict(2 => n -1)
        if optimize
            for solver in solvers
                @testset let n=n, solver=solver
                    @test poly_optimize(solver, sp).objective ≈ 1 atol = 1e-4 skip=solver==:SCSMoment
                end
            end
        end
    end
end

@testset "Example 6.1 Generalized Rosebrock function from correlative sparsity paper" begin
    for n in 600:100:1000
        DynamicPolynomials.@polyvar x[1:n]
        sp = Relaxation.SparsityCorrelative(
            poly_problem(1 + sum(100((x[i] - x[i-1]^2)^2 + (1 - x[i])^2) for i in 2:n))
        )
        @test StatsBase.countmap(length.(groupings(sp).var_cliques)) == Dict(2 => n -1)
        if optimize
            for solver in solvers
                @testset let n=n, solver=solver
                    @test(poly_optimize(solver, sp).objective ≈ 1, atol=solver==:HypatiaMoment ? 5e-2 : 3e-4,
                        skip=solver ∈ (:ClarabelMoment, :SCSMoment))
                    # Clarabel very inaccurate, SCS very slow
                end
            end
        end
    end
end

@testset "Example 6.2 Problem (6.1) from correlative sparsity paper" begin
    # Note that these values are generated with Mathematica. They do not correspond to the values given in the paper; however,
    # there, they use a chordal extension heuristic (an option that we also have, but it is a different one). Here, we disable
    # chordal completion in order to obtain the best possible (and unique) answer.
    # The numerical results for μ = 0 are exact, as NMinimize gives the same value for feasible points.
    results = Dict(
        (2, 4, 0., 6) => (Dict(3 => 1, 4 => 15, 5 => 15, 6 => 7), 0.036193),
        (2, 4, 0., 12) => (Dict(3 => 1, 4 => 39, 5 => 39, 6 => 19), 0.0440061),
        (2, 4, 0., 18) => (Dict(3 => 1, 4 => 63, 5 => 63, 6 => 31), 0.0517071),
        (2, 4, 0., 24) => (Dict(3 => 1, 4 => 87, 5 => 87, 6 => 43), 0.0594),
        (2, 4, 0., 30) => (Dict(3 => 1, 4 => 111, 5 => 111, 6 => 55), 0.0670932), # here Mosek has trouble when presolve is on
        (2, 4, .5, 6) => (Dict(4 => 3, 7 => 4, 10 => 3), 0.059118),
        (2, 4, .5, 12) => (Dict(4 => 3, 7 => 4, 10 => 9), 0.112406), # gap to NMinimize: 0.000079
        (2, 4, .5, 18) => (Dict(4 => 3, 7 => 4, 10 => 15), 0.165774), # gap to NMinimize: 0.00024
        (2, 4, .5, 24) => (Dict(4 => 3, 7 => 4, 10 => 21), 0.219139), # gap to NMinimize: 0.0004
        (2, 4, .5, 30) => (Dict(4 => 3, 7 => 4, 10 => 27), 0.272505), # gap to NMinimize: 0.00056
        (1, 2, 1, 6) => (Dict(2 => 1, 4 => 2, 5 => 3), 0.0282973),
        (1, 2, 1, 12) => (Dict(2 => 1, 4 => 2, 5 => 9), 0.052990), # gap to NMinimize: 0.0000019
        (1, 2, 1, 18) => (Dict(2 => 1, 4 => 2, 5 => 15), 0.077680), # gap to NMinimize: 0.0000054
        (1, 2, 1, 24) => (Dict(2 => 1, 4 => 2, 5 => 21), 0.102370), # gap to NMinimize: 0.000009
        (1, 2, 1, 30) => (Dict(2 => 1, 4 => 2, 5 => 27), 0.127060), # gap to NMinimize: 0.000012
    )
    for (nx, ny, μ) in ((2, 4, 0.), (2, 4, 0.5), (1, 2, 1))
        for M in 6:6:30
            # they don't add y₁ = 0 but instead make this implicit
            DynamicPolynomials.@polyvar x[1:M-1, 1:nx] y[2:M, 1:ny]
            ninv = inv(ny + nx)
            Y(i, j) = isone(i) ? 0 : y[i-1, j]
            sp = Relaxation.SparsityCorrelative(
                poly_problem(
                    sum(sum((Y(i, j) + .25)^4 for j in 1:ny) + sum((x[i, j] + .25)^4 for j in 1:nx)
                        for i in 1:M-1) +
                    sum((Y(M, j) + .25)^4 for j in 1:ny),
                    zero=[
                        sum((k == j ? 0.5 : (k == j +1 ? .25 : (k == j -1 ? -.25 : 0.))) #= A[j, k] =# * Y(i, k) for k in 1:ny) +
                        sum((j - k) * ninv #= B[j, k] =# * x[i, k] for k in 1:nx) +
                        sum(Y(i, k) * μ * (k + l) * ninv #= C[k, l] =# * x[i, l] for k in 1:ny for l in 1:nx) -
                        Y(i +1, j)
                        for i in 1:M-1 for j in 1:ny
                    ]
                ),
                chordal_completion=false
            )
            cliques, bound = results[(nx, ny, μ, M)]
            @test StatsBase.countmap(length.(groupings(sp).var_cliques)) == cliques
            if optimize && μ != 0.5
                for solver in solvers # the largest 0.5 can take about 180s to solve
                    if solver == :MosekSOS
                        parameters = Dict(:MSK_IPAR_PRESOLVE_USE => Mosek.MSK_PRESOLVE_MODE_OFF.value)
                    else
                        parameters = Dict()
                    end
                    @testset let nx=nx, ny=ny, μ=μ, M=M, solver=solver
                        @test poly_optimize(solver, sp; parameters...).objective ≈ bound atol = solver===:SCSMoment ? 1e-3 : 1e-4
                    end
                end
            end
        end
    end
end

@testset "Example 6.2 Problem (6.2) from correlative sparsity paper" begin
    for (M, result) in ((600, 0.00645), (700, 0.00553), (800, 0.00484), (900, 0.0043), (1000, 0.0038))
        DynamicPolynomials.@polyvar x[1:M] y[2:M]
        Y(i) = isone(i) ? 1 : y[i-1]
        sp = Relaxation.SparsityCorrelative(
            poly_problem(
                sum(Y(i)^2 + x[i]^2 for i in 1:M-1) / M,
                zero=[Y(i) + (Y(i)^2 - x[i]) - Y(i+1) for i in 1:M-1]
            )
        )
        # We use a more sophisticated analysis that the correlative sparsity paper: relaxation order 1 means that the
        # constraints will only allow for the constant prefactor - therefore, we can also count variables occurring
        # simultaneously in terms instead of the whole constraints!
        @test StatsBase.countmap(length.(groupings(sp).var_cliques)) == Dict(1 => 2*(M -1))
        if optimize
            for solver in solvers
                @testset let M=M, solver=solver
                    @test poly_optimize(solver, sp).objective ≈ result atol = 1e-4
                end
            end
        end
    end
end

@testset "Example 3.1 from correlative term sparsity paper" begin
    DynamicPolynomials.@polyvar x[1:3]
    sp = Relaxation.SparsityCorrelative(poly_problem(1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] * x[2] + x[2] * x[3] + x[3]), 2)
    @test strRep(sp) == "Relaxation.SparsityCorrelative of a polynomial optimization problem
Variable cliques:
  x[1], x[2]
  x[2], x[3]
PSD block sizes:
  [6 => 2]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.625 atol = solver===:SCSMoment ? 1e-5 : 1e-6
            end
        end
    end
end

@testset "Example 3.4 from correlative term sparsity paper" begin
    DynamicPolynomials.@polyvar x[1:6]
    sp = Relaxation.SparsityCorrelative(
        poly_problem(1 + sum(x[i]^4 for i in 1:6) + x[1] * x[2] * x[3] + x[3] * x[4] * x[5] + x[3] * x[4] * x[6] +
                     x[3] * x[5] * x[6] + x[4] * x[5] * x[6], perturbation=1e-4), 2
    )
    @test strRep(sp) == "Relaxation.SparsityCorrelative of a polynomial optimization problem
Variable cliques:
  x[3], x[4], x[5], x[6]
  x[1], x[2], x[3]
PSD block sizes:
  [15 => 1, 10 => 1]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.5042 atol = 1e-3
            end
        end
    end
end

@testset "Example 4.1 from Zhen, Fantuzzi, Papachristodoulou review" begin
    DynamicPolynomials.@polyvar x[1:50]
    sp = Relaxation.SparsityCorrelative(
        poly_problem(sum((x[i-1] + x[i] + x[i+1])^4 for i in 2:49)), 2
    )
    @test strRep(sp) == "Relaxation.SparsityCorrelative of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[3]
  x[2], x[3], x[4]
  x[3], x[4], x[5]
  x[4], x[5], x[6]
  x[5], x[6], x[7]
  x[6], x[7], x[8]
  x[7], x[8], x[9]
  x[8], x[9], x[10]
  x[9], x[10], x[11]
  x[10], x[11], x[12]
  x[11], x[12], x[13]
  x[12], x[13], x[14]
  x[13], x[14], x[15]
  x[14], x[15], x[16]
  x[15], x[16], x[17]
  x[16], x[17], x[18]
  x[17], x[18], x[19]
  x[18], x[19], x[20]
  x[19], x[20], x[21]
  x[20], x[21], x[22]
  x[21], x[22], x[23]
  x[22], x[23], x[24]
  x[23], x[24], x[25]
  x[24], x[25], x[26]
  x[25], x[26], x[27]
  x[26], x[27], x[28]
  x[27], x[28], x[29]
  x[28], x[29], x[30]
  x[29], x[30], x[31]
  x[30], x[31], x[32]
  x[31], x[32], x[33]
  x[32], x[33], x[34]
  x[33], x[34], x[35]
  x[34], x[35], x[36]
  x[35], x[36], x[37]
  x[36], x[37], x[38]
  x[37], x[38], x[39]
  x[38], x[39], x[40]
  x[39], x[40], x[41]
  x[40], x[41], x[42]
  x[41], x[42], x[43]
  x[42], x[43], x[44]
  x[43], x[44], x[45]
  x[44], x[45], x[46]
  x[45], x[46], x[47]
  x[46], x[47], x[48]
  x[47], x[48], x[49]
  x[48], x[49], x[50]
PSD block sizes:
  [10 => 48]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0 atol = solver===:SCSMoment ? 1e-5 : 1e-6
            end
        end
    end
end

@testset "Example 4.4 from Zhen, Fantuzzi, Papachristodoulou review" begin
    DynamicPolynomials.@polyvar x[1:4]
    prob = poly_problem(2 + x[1]^2 * x[4]^2 * (x[1]^2 * x[4]^2 - 1) - x[1]^2 + x[1]^4 +
                        sum((x[i]^2 * x[i-1]^2 * (x[i]^2 * x[i-1]^2 - 1) - x[i]^2 + x[i]^4 for i in 2:4)),
                        perturbation=1e-4)
    sp = Relaxation.SparsityCorrelative(prob, 4)
    @test strRep(sp) == "Relaxation.SparsityCorrelative of a polynomial optimization problem
Variable cliques:
  x[1], x[2], x[4]
  x[2], x[3], x[4]
PSD block sizes:
  [35 => 2]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.11008 atol = 1e-3
            end
        end
    end

    sp = Relaxation.SparsityCorrelative(prob, 4, chordal_completion=false)
    @test strRep(sp) == "Relaxation.SparsityCorrelative of a polynomial optimization problem
Variable cliques:
  x[1], x[2]
  x[1], x[4]
  x[2], x[3]
  x[3], x[4]
PSD block sizes:
  [15 => 4]"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, sp).objective ≈ 0.11008 atol = 1e-3
            end
        end
    end

    @test isnothing(iterate!(sp))
end

@testset "Example 6.1 from Josz, Molzahn" begin
    DynamicPolynomials.@polyvar x[1:4]
    prob = poly_problem(x[1]*x[2] + x[1]*x[4], nonneg=[x[1]*x[2]+x[1]*x[3], x[1]*x[3]+x[1]*x[4]+x[1]*x[2]])
    @test strRep(groupings(Relaxation.SparsityCorrelative(prob, 2, low_order_nonneg=[2], chordal_completion=false))) ==
        "Groupings for the relaxation of a polynomial optimization problem
Variable cliques
================
[x₁, x₂, x₃]
[x₁, x₄]

Block groupings
===============
Objective: 2 blocks
  10 [1, x₃, x₂, x₁, x₃², x₂x₃, x₂², x₁x₃, x₁x₂, x₁²]
   6 [1, x₄, x₁, x₄², x₁x₄, x₁²]
Nonnegative constraint #1: 1 block
  4 [1, x₃, x₂, x₁]
Nonnegative constraint #2: 1 block
  1 [1]"
end

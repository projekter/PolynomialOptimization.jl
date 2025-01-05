# These tests require Mosek!
include("./shared.jl")
using PolynomialOptimization: tightening_methods

function roundproblem(args...; kwargs...)
    roundfn = x -> round(x, digits=5)
    p = poly_problem(args...; kwargs...)
    map_coefficients!.(roundfn, p.objective)
    map_coefficients!.(roundfn, p.constr_zero)
    map_coefficients!.(roundfn, p.constr_nonneg)
    map_coefficients!.(roundfn, p.constr_psd)
    return p
end

@testset "Example 6.1" begin
    DynamicPolynomials.@polyvar x[1:3]
    local prob
    for method in tightening_methods
        prob = roundproblem(x[1] * x[2] * (10 - x[3]), nonneg=[x..., 1 - sum(x)], tighter=method)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: 10.0x₁x₂ - x₁x₂x₃
4 equality constraints
01: 10.0x₁x₂ - x₁x₂x₃ - 20.0x₁²x₂ + 3.0x₁²x₂x₃ = 0
02: 10.0x₁x₂ - x₁x₂x₃ - 20.0x₁x₂² + 3.0x₁x₂²x₃ = 0
03: -21.0x₁x₂x₃ + 3.0x₁x₂x₃² = 0
04: -20.0x₁x₂ + 23.0x₁x₂x₃ + 20.0x₁x₂² + 20.0x₁²x₂ - 3.0x₁x₂x₃² - 3.0x₁x₂²x₃ - 3.0x₁²x₂x₃ = 0
8 nonnegative constraints
05: x₁ ≥ 0
06: x₂ ≥ 0
07: x₃ ≥ 0
08: 1.0 - x₃ - x₂ - x₁ ≥ 0
09: 10.0x₂ - x₂x₃ - 20.0x₁x₂ + 3.0x₁x₂x₃ ≥ 0
10: 10.0x₁ - x₁x₃ - 20.0x₁x₂ + 3.0x₁x₂x₃ ≥ 0
11: -21.0x₁x₂ + 3.0x₁x₂x₃ ≥ 0
12: -20.0x₁x₂ + 3.0x₁x₂x₃ ≥ 0"
    end
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ -0.05208 atol = solver==:SCSMoment ? 1e-3 : 1e-5
            end
        end
    end
end

@testset "Example 6.2" begin
    DynamicPolynomials.@polyvar x[1:3]
    local prob
    for method in tightening_methods
        prob = roundproblem(x[1]^4 * x[2]^2 + x[1]^2 * x[2]^4 + x[3]^6 - 3prod(x .^ 2) + sum(x .^ 4), nonneg=[sum(x .^ 2)-1],
            tighter=method)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: x₃⁴ + x₂⁴ + x₁⁴ + x₃⁶ - 3.0x₁²x₂²x₃² + x₁²x₂⁴ + x₁⁴x₂²
4 equality constraints
1: 4.0x₁³ - 4.0x₁x₃⁴ - 6.0x₁x₂²x₃² - 2.0x₁x₂⁴ + 4.0x₁³x₂² - 4.0x₁⁵ - 6.0x₁x₃⁶ + 18.0x₁³x₂²x₃² - 6.0x₁³x₂⁴ - 6.0x₁⁵x₂² = 0
2: 4.0x₂³ - 4.0x₂x₃⁴ - 4.0x₂⁵ - 6.0x₁²x₂x₃² + 4.0x₁²x₂³ - 2.0x₁⁴x₂ - 6.0x₂x₃⁶ + 18.0x₁²x₂³x₃² - 6.0x₁²x₂⁵ - 6.0x₁⁴x₂³ = 0
3: 4.0x₃³ + 2.0x₃⁵ - 4.0x₂⁴x₃ - 6.0x₁²x₂²x₃ - 4.0x₁⁴x₃ - 6.0x₃⁷ + 18.0x₁²x₂²x₃³ - 6.0x₁²x₂⁴x₃ - 6.0x₁⁴x₂²x₃ = 0
4: -2.0x₃⁴ - 2.0x₂⁴ - 2.0x₁⁴ - x₃⁶ + 2.0x₂²x₃⁴ + 2.0x₂⁴x₃² + 2.0x₂⁶ + 2.0x₁²x₃⁴ + 9.0x₁²x₂²x₃² - x₁²x₂⁴ + 2.0x₁⁴x₃² - x₁⁴x₂² + 2.0x₁⁶ + 3.0x₃⁸ + 3.0x₂²x₃⁶ + 3.0x₁²x₃⁶ - 9.0x₁²x₂²x₃⁴ - 6.0x₁²x₂⁴x₃² + 3.0x₁²x₂⁶ - 6.0x₁⁴x₂²x₃² + 6.0x₁⁴x₂⁴ + 3.0x₁⁶x₂² = 0
2 nonnegative constraints
5: -1.0 + x₃² + x₂² + x₁² ≥ 0
6: 2.0x₃⁴ + 2.0x₂⁴ + 2.0x₁⁴ + 3.0x₃⁶ - 9.0x₁²x₂²x₃² + 3.0x₁²x₂⁴ + 3.0x₁⁴x₂² ≥ 0"
    end
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 1/3 atol = solver==:SCSMoment ? 1e-3 : 1e-6
            end
        end
    end
end

@testset "Example 6.3" begin
    DynamicPolynomials.@polyvar x[1:4]
    local prob
    for method in tightening_methods
        prob = roundproblem(x[1]*x[2] + x[2]*x[3] + x[3]*x[4] - 3prod(x) + sum(x .^ 3),
            nonneg=[x..., 1 - x[1] - x[2], 1 - x[3] - x[4]], tighter=method)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 4 variables
Objective: x₃x₄ + x₂x₃ + x₁x₂ + x₄³ + x₃³ + x₂³ + x₁³ - 3.0x₁x₂x₃x₄
6 equality constraints
01: x₁x₂ - x₁x₂x₃ - 2.0x₁²x₂ + 3.0x₁³ - 3.0x₁x₂x₃x₄ - 3.0x₁x₂³ - 3.0x₁⁴ + 6.0x₁²x₂x₃x₄ = 0
02: x₂x₃ + x₁x₂ - x₂²x₃ + 3.0x₂³ - 2.0x₁x₂² - 3.0x₂⁴ - 3.0x₁x₂x₃x₄ - 3.0x₁³x₂ + 6.0x₁x₂²x₃x₄ = 0
03: x₃x₄ + x₂x₃ - 2.0x₃²x₄ + 3.0x₃³ - x₂x₃² - 3.0x₃x₄³ - 3.0x₃⁴ - 3.0x₁x₂x₃x₄ + 6.0x₁x₂x₃²x₄ = 0
04: x₃x₄ + 3.0x₄³ - 2.0x₃x₄² - x₂x₃x₄ - 3.0x₄⁴ - 3.0x₃³x₄ - 3.0x₁x₂x₃x₄ + 6.0x₁x₂x₃x₄² = 0
05: -x₂x₃ - 2.0x₁x₂ + x₂²x₃ - 3.0x₂³ + x₁x₂x₃ + 2.0x₁x₂² + 2.0x₁²x₂ - 3.0x₁³ + 3.0x₂⁴ + 6.0x₁x₂x₃x₄ + 3.0x₁x₂³ + 3.0x₁³x₂ + 3.0x₁⁴ - 6.0x₁x₂²x₃x₄ - 6.0x₁²x₂x₃x₄ = 0
06: -2.0x₃x₄ - x₂x₃ - 3.0x₄³ + 2.0x₃x₄² + 2.0x₃²x₄ - 3.0x₃³ + x₂x₃x₄ + x₂x₃² + 3.0x₄⁴ + 3.0x₃x₄³ + 3.0x₃³x₄ + 3.0x₃⁴ + 6.0x₁x₂x₃x₄ - 6.0x₁x₂x₃x₄² - 6.0x₁x₂x₃²x₄ = 0
12 nonnegative constraints
07: x₁ ≥ 0
08: x₂ ≥ 0
09: x₃ ≥ 0
10: x₄ ≥ 0
11: 1.0 - x₂ - x₁ ≥ 0
12: 1.0 - x₄ - x₃ ≥ 0
13: x₂ - x₂x₃ - 2.0x₁x₂ + 3.0x₁² - 3.0x₂x₃x₄ - 3.0x₂³ - 3.0x₁³ + 6.0x₁x₂x₃x₄ ≥ 0
14: x₃ + x₁ - x₂x₃ + 3.0x₂² - 2.0x₁x₂ - 3.0x₂³ - 3.0x₁x₃x₄ - 3.0x₁³ + 6.0x₁x₂x₃x₄ ≥ 0
15: x₄ + x₂ - 2.0x₃x₄ + 3.0x₃² - x₂x₃ - 3.0x₄³ - 3.0x₃³ - 3.0x₁x₂x₄ + 6.0x₁x₂x₃x₄ ≥ 0
16: x₃ + 3.0x₄² - 2.0x₃x₄ - x₂x₃ - 3.0x₄³ - 3.0x₃³ - 3.0x₁x₂x₃ + 6.0x₁x₂x₃x₄ ≥ 0
17: -x₂x₃ - 2.0x₁x₂ - 3.0x₂³ - 3.0x₁³ + 6.0x₁x₂x₃x₄ ≥ 0
18: -2.0x₃x₄ - x₂x₃ - 3.0x₄³ - 3.0x₃³ + 6.0x₁x₂x₃x₄ ≥ 0"
    end
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = solver==:SCSMoment ? 1e-3 : 1e-5
            end
        end
    end
end

@testset "Example 6.4" begin
    DynamicPolynomials.@polyvar x[1:2]
    local prob
    for method in tightening_methods
        prob = roundproblem(x[1]^2 + 50.0x[2]^2, nonneg=[x[1]^2-.5, x[2]^2-2x[1]*x[2]-.125, x[2]^2+2x[1]*x[2]-0.125],
            tighter=method)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 50.0x₂² + x₁²
5 equality constraints
01: 2.0x₁ + 18.4x₁x₂² - 0.8x₁³ - 26035.2x₁x₂⁴ + 323.2x₁³x₂² - 6.4x₁⁵ + 11520.0x₁x₂⁶ + 230.4x₁³x₂⁴ = 0
02: 100.0x₂ - 880.0x₂³ - 14.4x₁²x₂ + 640.0x₂⁵ - 27155.2x₁²x₂³ + 28.8x₁⁴x₂ + 11520.0x₁²x₂⁵ + 230.4x₁⁴x₂³ = 0
03: 5.0x₂² - 0.2x₁² - 40.0x₂⁴ + 70.0x₁²x₂² - 1.2x₁⁴ + 80.0x₁²x₂⁴ - 160.0x₁⁴x₂² + 3.2x₁⁶ = 0
04: -27.5x₂² + 0.05x₁x₂ - 0.4x₁² + 240.0x₂⁴ + 368.2x₁x₂³ - 36.0x₁²x₂² - 6.5x₁³x₂ + 0.8x₁⁴ - 160.0x₂⁶ - 6508.8x₁x₂⁵ + 13257.6x₁²x₂⁴ - 646.4x₁³x₂³ - 8.0x₁⁴x₂² + 12.8x₁⁵x₂ + 2880.0x₁x₂⁷ - 5760.0x₁²x₂⁶ + 57.6x₁³x₂⁵ - 115.2x₁⁴x₂⁴ = 0
05: -27.5x₂² - 0.05x₁x₂ - 0.4x₁² + 240.0x₂⁴ - 368.2x₁x₂³ - 36.0x₁²x₂² + 6.5x₁³x₂ + 0.8x₁⁴ - 160.0x₂⁶ + 6508.8x₁x₂⁵ + 13257.6x₁²x₂⁴ + 646.4x₁³x₂³ - 8.0x₁⁴x₂² - 12.8x₁⁵x₂ - 2880.0x₁x₂⁷ - 5760.0x₁²x₂⁶ - 57.6x₁³x₂⁵ - 115.2x₁⁴x₂⁴ = 0
6 nonnegative constraints
06: -0.5 + x₁² ≥ 0
07: -0.125 + x₂² - 2.0x₁x₂ ≥ 0
08: -0.125 + x₂² + 2.0x₁x₂ ≥ 0
09: -10.0x₂² + 0.4x₁² + 80.0x₂⁴ - 160.0x₁²x₂² + 3.2x₁⁴ ≥ 0
10: 220.0x₂² - 0.4x₁x₂ + 3.2x₁² - 160.0x₂⁴ - 6468.8x₁x₂³ + 320.0x₁²x₂² + 0.8x₁³x₂ - 6.4x₁⁴ + 2880.0x₁x₂⁵ + 57.6x₁³x₂³ ≥ 0
11: 220.0x₂² + 0.4x₁x₂ + 3.2x₁² - 160.0x₂⁴ + 6468.8x₁x₂³ + 320.0x₁²x₂² - 0.8x₁³x₂ - 6.4x₁⁴ - 2880.0x₁x₂⁵ - 57.6x₁³x₂³ ≥ 0"
    end
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 112.6516994 atol = 1e-2
            end
        end
    end
end

@testset "Example 6.5" begin
    DynamicPolynomials.@polyvar x[1:3]
    local prob
    for method in tightening_methods
        prob = roundproblem(sum(x .^ 3) + 4prod(x) - (x[1]*(x[2]^2+x[3]^2) + x[2]*(x[3]^2+x[1]^2) + x[3]*(x[1]^2+x[2]^2)),
            nonneg=[x[1], x[1]*x[2]-1, x[2]*x[3]-1], tighter=method)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: x₃³ - x₂x₃² - x₂²x₃ + x₂³ - x₁x₃² + 4.0x₁x₂x₃ - x₁x₂² - x₁²x₃ - x₁²x₂ + x₁³
5 equality constraints
01: -x₃² - 2.0x₂x₃ + 3.0x₂² + 4.0x₁x₃ - 2.0x₁x₂ - x₁² - 3.0x₃⁴ + 2.0x₂x₃³ + x₂²x₃² + 5.0x₁x₃³ - 5.0x₁x₂x₃² + x₁x₂²x₃ - 3.0x₁x₂³ - x₁²x₃² + 2.0x₁²x₂² - x₁³x₃ + x₁³x₂ = 0
02: 3.0x₃² - 2.0x₂x₃ - x₂² - 2.0x₁x₃ + 4.0x₁x₂ - x₁² - 3.0x₂x₃³ + 2.0x₂²x₃² + x₂³x₃ + 2.0x₁x₂x₃² - 4.0x₁x₂²x₃ + x₁²x₂x₃ = 0
03: -x₁x₃² + 4.0x₁x₂x₃ - x₁x₂² - 2.0x₁²x₃ - 2.0x₁²x₂ + 3.0x₁³ + 3.0x₁x₂x₃³ - x₁x₂²x₃² + x₁x₂³x₃ - 3.0x₁x₂⁴ - 2.0x₁²x₂x₃² + 2.0x₁²x₂³ - x₁³x₂x₃ + x₁³x₂² = 0
04: 3.0x₃³ - x₂x₃² + x₂²x₃ - 3.0x₂³ - 2.0x₁x₃² + 2.0x₁x₂² - x₁²x₃ + x₁²x₂ - 3.0x₁x₂x₃³ + x₁x₂²x₃² - x₁x₂³x₃ + 3.0x₁x₂⁴ + 2.0x₁²x₂x₃² - 2.0x₁²x₂³ + x₁³x₂x₃ - x₁³x₂² = 0
05: -3.0x₃³ + 2.0x₂x₃² + x₂²x₃ + 2.0x₁x₃² - 4.0x₁x₂x₃ + x₁²x₃ + 3.0x₂x₃⁴ - 2.0x₂²x₃³ - x₂³x₃² - 2.0x₁x₂x₃³ + 4.0x₁x₂²x₃² - x₁²x₂x₃² = 0
6 nonnegative constraints
06: x₁ ≥ 0
07: -1.0 + x₁x₂ ≥ 0
08: -1.0 + x₂x₃ ≥ 0
09: -x₃² + 4.0x₂x₃ - x₂² - 2.0x₁x₃ - 2.0x₁x₂ + 3.0x₁² + 3.0x₂x₃³ - x₂²x₃² + x₂³x₃ - 3.0x₂⁴ - 2.0x₁x₂x₃² + 2.0x₁x₂³ - x₁²x₂x₃ + x₁²x₂² ≥ 0
10: -3.0x₃³ + x₂x₃² - x₂²x₃ + 3.0x₂³ + 2.0x₁x₃² - 2.0x₁x₂² + x₁²x₃ - x₁²x₂ ≥ 0
11: 3.0x₃³ - 2.0x₂x₃² - x₂²x₃ - 2.0x₁x₃² + 4.0x₁x₂x₃ - x₁²x₃ ≥ 0"
    end
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob, 4).objective ≈ 0.9491545329 atol = solver==:SCSMoment ? 1e-4 : 1e-6
                # deg 3 is possible, but unstable
            end
        end
    end
end

@testset "Example 6.6" begin
    DynamicPolynomials.@polyvar x[1:4]
    X(i) = i == 0 ? 1 : x[i]
    local prob
    for method in tightening_methods
        prob = roundproblem(sum(x .^ 2) + sum(prod(i == j ? 1 : X(i) - X(j) for j in 0:4) for i in 0:4), nonneg=x.^2 .- 1,
            tighter=method)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 4 variables
Objective: 1.0 - x₄ - x₃ - x₂ - x₁ + x₄² + x₃x₄ + x₃² + x₂x₄ + x₂x₃ + x₂² + x₁x₄ + x₁x₃ + x₁x₂ + x₁² - x₄³ + x₃x₄² + x₃²x₄ - x₃³ + x₂x₄² - 3.0x₂x₃x₄ + x₂x₃² + x₂²x₄ + x₂²x₃ - x₂³ + x₁x₄² - 3.0x₁x₃x₄ + x₁x₃² - 3.0x₁x₂x₄ - 3.0x₁x₂x₃ + x₁x₂² + x₁²x₄ + x₁²x₃ + x₁²x₂ - x₁³ + x₄⁴ - x₃x₄³ - x₃³x₄ + x₃⁴ - x₂x₄³ + x₂x₃x₄² + x₂x₃²x₄ - x₂x₃³ + x₂²x₃x₄ - x₂³x₄ - x₂³x₃ + x₂⁴ - x₁x₄³ + x₁x₃x₄² + x₁x₃²x₄ - x₁x₃³ + x₁x₂x₄² - 3.0x₁x₂x₃x₄ + x₁x₂x₃² + x₁x₂²x₄ + x₁x₂²x₃ - x₁x₂³ + x₁²x₃x₄ + x₁²x₂x₄ + x₁²x₂x₃ - x₁³x₄ - x₁³x₃ - x₁³x₂ + x₁⁴
8 equality constraints
01: -1.0 + x₄ + x₃ + x₂ + 2.0x₁ + x₄² - 3.0x₃x₄ + x₃² - 3.0x₂x₄ - 3.0x₂x₃ + x₂² + 2.0x₁x₄ + 2.0x₁x₃ + 2.0x₁x₂ - 2.0x₁² - x₄³ + x₃x₄² + x₃²x₄ - x₃³ + x₂x₄² - 3.0x₂x₃x₄ + x₂x₃² + x₂²x₄ + x₂²x₃ - x₂³ + 2.0x₁x₃x₄ + 2.0x₁x₂x₄ + 2.0x₁x₂x₃ - 4.0x₁²x₄ - 4.0x₁²x₃ - 4.0x₁²x₂ + 2.0x₁³ - x₁²x₄² + 3.0x₁²x₃x₄ - x₁²x₃² + 3.0x₁²x₂x₄ + 3.0x₁²x₂x₃ - x₁²x₂² - 2.0x₁³x₄ - 2.0x₁³x₃ - 2.0x₁³x₂ + 3.0x₁⁴ + x₁²x₄³ - x₁²x₃x₄² - x₁²x₃²x₄ + x₁²x₃³ - x₁²x₂x₄² + 3.0x₁²x₂x₃x₄ - x₁²x₂x₃² - x₁²x₂²x₄ - x₁²x₂²x₃ + x₁²x₂³ - 2.0x₁³x₃x₄ - 2.0x₁³x₂x₄ - 2.0x₁³x₂x₃ + 3.0x₁⁴x₄ + 3.0x₁⁴x₃ + 3.0x₁⁴x₂ - 4.0x₁⁵ = 0
02: -1.0 + x₄ + x₃ + 2.0x₂ + x₁ + x₄² - 3.0x₃x₄ + x₃² + 2.0x₂x₄ + 2.0x₂x₃ - 2.0x₂² - 3.0x₁x₄ - 3.0x₁x₃ + 2.0x₁x₂ + x₁² - x₄³ + x₃x₄² + x₃²x₄ - x₃³ + 2.0x₂x₃x₄ - 4.0x₂²x₄ - 4.0x₂²x₃ + 2.0x₂³ + x₁x₄² - 3.0x₁x₃x₄ + x₁x₃² + 2.0x₁x₂x₄ + 2.0x₁x₂x₃ - 4.0x₁x₂² + x₁²x₄ + x₁²x₃ - x₁³ - x₂²x₄² + 3.0x₂²x₃x₄ - x₂²x₃² - 2.0x₂³x₄ - 2.0x₂³x₃ + 3.0x₂⁴ + 3.0x₁x₂²x₄ + 3.0x₁x₂²x₃ - 2.0x₁x₂³ - x₁²x₂² + x₂²x₄³ - x₂²x₃x₄² - x₂²x₃²x₄ + x₂²x₃³ - 2.0x₂³x₃x₄ + 3.0x₂⁴x₄ + 3.0x₂⁴x₃ - 4.0x₂⁵ - x₁x₂²x₄² + 3.0x₁x₂²x₃x₄ - x₁x₂²x₃² - 2.0x₁x₂³x₄ - 2.0x₁x₂³x₃ + 3.0x₁x₂⁴ - x₁²x₂²x₄ - x₁²x₂²x₃ + x₁³x₂² = 0
03: -1.0 + x₄ + 2.0x₃ + x₂ + x₁ + x₄² + 2.0x₃x₄ - 2.0x₃² - 3.0x₂x₄ + 2.0x₂x₃ + x₂² - 3.0x₁x₄ + 2.0x₁x₃ - 3.0x₁x₂ + x₁² - x₄³ - 4.0x₃²x₄ + 2.0x₃³ + x₂x₄² + 2.0x₂x₃x₄ - 4.0x₂x₃² + x₂²x₄ - x₂³ + x₁x₄² + 2.0x₁x₃x₄ - 4.0x₁x₃² - 3.0x₁x₂x₄ + 2.0x₁x₂x₃ + x₁x₂² + x₁²x₄ + x₁²x₂ - x₁³ - x₃²x₄² - 2.0x₃³x₄ + 3.0x₃⁴ + 3.0x₂x₃²x₄ - 2.0x₂x₃³ - x₂²x₃² + 3.0x₁x₃²x₄ - 2.0x₁x₃³ + 3.0x₁x₂x₃² - x₁²x₃² + x₃²x₄³ + 3.0x₃⁴x₄ - 4.0x₃⁵ - x₂x₃²x₄² - 2.0x₂x₃³x₄ + 3.0x₂x₃⁴ - x₂²x₃²x₄ + x₂³x₃² - x₁x₃²x₄² - 2.0x₁x₃³x₄ + 3.0x₁x₃⁴ + 3.0x₁x₂x₃²x₄ - 2.0x₁x₂x₃³ - x₁x₂²x₃² - x₁²x₃²x₄ - x₁²x₂x₃² + x₁³x₃² = 0
04: -1.0 + 2.0x₄ + x₃ + x₂ + x₁ - 2.0x₄² + 2.0x₃x₄ + x₃² + 2.0x₂x₄ - 3.0x₂x₃ + x₂² + 2.0x₁x₄ - 3.0x₁x₃ - 3.0x₁x₂ + x₁² + 2.0x₄³ - 4.0x₃x₄² - x₃³ - 4.0x₂x₄² + 2.0x₂x₃x₄ + x₂x₃² + x₂²x₃ - x₂³ - 4.0x₁x₄² + 2.0x₁x₃x₄ + x₁x₃² + 2.0x₁x₂x₄ - 3.0x₁x₂x₃ + x₁x₂² + x₁²x₃ + x₁²x₂ - x₁³ + 3.0x₄⁴ - 2.0x₃x₄³ - x₃²x₄² - 2.0x₂x₄³ + 3.0x₂x₃x₄² - x₂²x₄² - 2.0x₁x₄³ + 3.0x₁x₃x₄² + 3.0x₁x₂x₄² - x₁²x₄² - 4.0x₄⁵ + 3.0x₃x₄⁴ + x₃³x₄² + 3.0x₂x₄⁴ - 2.0x₂x₃x₄³ - x₂x₃²x₄² - x₂²x₃x₄² + x₂³x₄² + 3.0x₁x₄⁴ - 2.0x₁x₃x₄³ - x₁x₃²x₄² - 2.0x₁x₂x₄³ + 3.0x₁x₂x₃x₄² - x₁x₂²x₄² - x₁²x₃x₄² - x₁²x₂x₄² + x₁³x₄² = 0
05: 0.5x₁ - 0.5x₁x₄ - 0.5x₁x₃ - 0.5x₁x₂ - x₁² - 0.5x₁x₄² + 1.5x₁x₃x₄ - 0.5x₁x₃² + 1.5x₁x₂x₄ + 1.5x₁x₂x₃ - 0.5x₁x₂² - x₁²x₄ - x₁²x₃ - x₁²x₂ + x₁³ + 0.5x₁x₄³ - 0.5x₁x₃x₄² - 0.5x₁x₃²x₄ + 0.5x₁x₃³ - 0.5x₁x₂x₄² + 1.5x₁x₂x₃x₄ - 0.5x₁x₂x₃² - 0.5x₁x₂²x₄ - 0.5x₁x₂²x₃ + 0.5x₁x₂³ - x₁²x₃x₄ - x₁²x₂x₄ - x₁²x₂x₃ + 2.0x₁³x₄ + 2.0x₁³x₃ + 2.0x₁³x₂ - x₁⁴ + 0.5x₁³x₄² - 1.5x₁³x₃x₄ + 0.5x₁³x₃² - 1.5x₁³x₂x₄ - 1.5x₁³x₂x₃ + 0.5x₁³x₂² + x₁⁴x₄ + x₁⁴x₃ + x₁⁴x₂ - 1.5x₁⁵ - 0.5x₁³x₄³ + 0.5x₁³x₃x₄² + 0.5x₁³x₃²x₄ - 0.5x₁³x₃³ + 0.5x₁³x₂x₄² - 1.5x₁³x₂x₃x₄ + 0.5x₁³x₂x₃² + 0.5x₁³x₂²x₄ + 0.5x₁³x₂²x₃ - 0.5x₁³x₂³ + x₁⁴x₃x₄ + x₁⁴x₂x₄ + x₁⁴x₂x₃ - 1.5x₁⁵x₄ - 1.5x₁⁵x₃ - 1.5x₁⁵x₂ + 2.0x₁⁶ = 0
06: 0.5x₂ - 0.5x₂x₄ - 0.5x₂x₃ - x₂² - 0.5x₁x₂ - 0.5x₂x₄² + 1.5x₂x₃x₄ - 0.5x₂x₃² - x₂²x₄ - x₂²x₃ + x₂³ + 1.5x₁x₂x₄ + 1.5x₁x₂x₃ - x₁x₂² - 0.5x₁²x₂ + 0.5x₂x₄³ - 0.5x₂x₃x₄² - 0.5x₂x₃²x₄ + 0.5x₂x₃³ - x₂²x₃x₄ + 2.0x₂³x₄ + 2.0x₂³x₃ - x₂⁴ - 0.5x₁x₂x₄² + 1.5x₁x₂x₃x₄ - 0.5x₁x₂x₃² - x₁x₂²x₄ - x₁x₂²x₃ + 2.0x₁x₂³ - 0.5x₁²x₂x₄ - 0.5x₁²x₂x₃ + 0.5x₁³x₂ + 0.5x₂³x₄² - 1.5x₂³x₃x₄ + 0.5x₂³x₃² + x₂⁴x₄ + x₂⁴x₃ - 1.5x₂⁵ - 1.5x₁x₂³x₄ - 1.5x₁x₂³x₃ + x₁x₂⁴ + 0.5x₁²x₂³ - 0.5x₂³x₄³ + 0.5x₂³x₃x₄² + 0.5x₂³x₃²x₄ - 0.5x₂³x₃³ + x₂⁴x₃x₄ - 1.5x₂⁵x₄ - 1.5x₂⁵x₃ + 2.0x₂⁶ + 0.5x₁x₂³x₄² - 1.5x₁x₂³x₃x₄ + 0.5x₁x₂³x₃² + x₁x₂⁴x₄ + x₁x₂⁴x₃ - 1.5x₁x₂⁵ + 0.5x₁²x₂³x₄ + 0.5x₁²x₂³x₃ - 0.5x₁³x₂³ = 0
07: 0.5x₃ - 0.5x₃x₄ - x₃² - 0.5x₂x₃ - 0.5x₁x₃ - 0.5x₃x₄² - x₃²x₄ + x₃³ + 1.5x₂x₃x₄ - x₂x₃² - 0.5x₂²x₃ + 1.5x₁x₃x₄ - x₁x₃² + 1.5x₁x₂x₃ - 0.5x₁²x₃ + 0.5x₃x₄³ + 2.0x₃³x₄ - x₃⁴ - 0.5x₂x₃x₄² - x₂x₃²x₄ + 2.0x₂x₃³ - 0.5x₂²x₃x₄ + 0.5x₂³x₃ - 0.5x₁x₃x₄² - x₁x₃²x₄ + 2.0x₁x₃³ + 1.5x₁x₂x₃x₄ - x₁x₂x₃² - 0.5x₁x₂²x₃ - 0.5x₁²x₃x₄ - 0.5x₁²x₂x₃ + 0.5x₁³x₃ + 0.5x₃³x₄² + x₃⁴x₄ - 1.5x₃⁵ - 1.5x₂x₃³x₄ + x₂x₃⁴ + 0.5x₂²x₃³ - 1.5x₁x₃³x₄ + x₁x₃⁴ - 1.5x₁x₂x₃³ + 0.5x₁²x₃³ - 0.5x₃³x₄³ - 1.5x₃⁵x₄ + 2.0x₃⁶ + 0.5x₂x₃³x₄² + x₂x₃⁴x₄ - 1.5x₂x₃⁵ + 0.5x₂²x₃³x₄ - 0.5x₂³x₃³ + 0.5x₁x₃³x₄² + x₁x₃⁴x₄ - 1.5x₁x₃⁵ - 1.5x₁x₂x₃³x₄ + x₁x₂x₃⁴ + 0.5x₁x₂²x₃³ + 0.5x₁²x₃³x₄ + 0.5x₁²x₂x₃³ - 0.5x₁³x₃³ = 0
08: 0.5x₄ - x₄² - 0.5x₃x₄ - 0.5x₂x₄ - 0.5x₁x₄ + x₄³ - x₃x₄² - 0.5x₃²x₄ - x₂x₄² + 1.5x₂x₃x₄ - 0.5x₂²x₄ - x₁x₄² + 1.5x₁x₃x₄ + 1.5x₁x₂x₄ - 0.5x₁²x₄ - x₄⁴ + 2.0x₃x₄³ + 0.5x₃³x₄ + 2.0x₂x₄³ - x₂x₃x₄² - 0.5x₂x₃²x₄ - 0.5x₂²x₃x₄ + 0.5x₂³x₄ + 2.0x₁x₄³ - x₁x₃x₄² - 0.5x₁x₃²x₄ - x₁x₂x₄² + 1.5x₁x₂x₃x₄ - 0.5x₁x₂²x₄ - 0.5x₁²x₃x₄ - 0.5x₁²x₂x₄ + 0.5x₁³x₄ - 1.5x₄⁵ + x₃x₄⁴ + 0.5x₃²x₄³ + x₂x₄⁴ - 1.5x₂x₃x₄³ + 0.5x₂²x₄³ + x₁x₄⁴ - 1.5x₁x₃x₄³ - 1.5x₁x₂x₄³ + 0.5x₁²x₄³ + 2.0x₄⁶ - 1.5x₃x₄⁵ - 0.5x₃³x₄³ - 1.5x₂x₄⁵ + x₂x₃x₄⁴ + 0.5x₂x₃²x₄³ + 0.5x₂²x₃x₄³ - 0.5x₂³x₄³ - 1.5x₁x₄⁵ + x₁x₃x₄⁴ + 0.5x₁x₃²x₄³ + x₁x₂x₄⁴ - 1.5x₁x₂x₃x₄³ + 0.5x₁x₂²x₄³ + 0.5x₁²x₃x₄³ + 0.5x₁²x₂x₄³ - 0.5x₁³x₄³ = 0
8 nonnegative constraints
09: -1.0 + x₁² ≥ 0
10: -1.0 + x₂² ≥ 0
11: -1.0 + x₃² ≥ 0
12: -1.0 + x₄² ≥ 0
13: -0.5x₁ + 0.5x₁x₄ + 0.5x₁x₃ + 0.5x₁x₂ + x₁² + 0.5x₁x₄² - 1.5x₁x₃x₄ + 0.5x₁x₃² - 1.5x₁x₂x₄ - 1.5x₁x₂x₃ + 0.5x₁x₂² + x₁²x₄ + x₁²x₃ + x₁²x₂ - 1.5x₁³ - 0.5x₁x₄³ + 0.5x₁x₃x₄² + 0.5x₁x₃²x₄ - 0.5x₁x₃³ + 0.5x₁x₂x₄² - 1.5x₁x₂x₃x₄ + 0.5x₁x₂x₃² + 0.5x₁x₂²x₄ + 0.5x₁x₂²x₃ - 0.5x₁x₂³ + x₁²x₃x₄ + x₁²x₂x₄ + x₁²x₂x₃ - 1.5x₁³x₄ - 1.5x₁³x₃ - 1.5x₁³x₂ + 2.0x₁⁴ ≥ 0
14: -0.5x₂ + 0.5x₂x₄ + 0.5x₂x₃ + x₂² + 0.5x₁x₂ + 0.5x₂x₄² - 1.5x₂x₃x₄ + 0.5x₂x₃² + x₂²x₄ + x₂²x₃ - 1.5x₂³ - 1.5x₁x₂x₄ - 1.5x₁x₂x₃ + x₁x₂² + 0.5x₁²x₂ - 0.5x₂x₄³ + 0.5x₂x₃x₄² + 0.5x₂x₃²x₄ - 0.5x₂x₃³ + x₂²x₃x₄ - 1.5x₂³x₄ - 1.5x₂³x₃ + 2.0x₂⁴ + 0.5x₁x₂x₄² - 1.5x₁x₂x₃x₄ + 0.5x₁x₂x₃² + x₁x₂²x₄ + x₁x₂²x₃ - 1.5x₁x₂³ + 0.5x₁²x₂x₄ + 0.5x₁²x₂x₃ - 0.5x₁³x₂ ≥ 0
15: -0.5x₃ + 0.5x₃x₄ + x₃² + 0.5x₂x₃ + 0.5x₁x₃ + 0.5x₃x₄² + x₃²x₄ - 1.5x₃³ - 1.5x₂x₃x₄ + x₂x₃² + 0.5x₂²x₃ - 1.5x₁x₃x₄ + x₁x₃² - 1.5x₁x₂x₃ + 0.5x₁²x₃ - 0.5x₃x₄³ - 1.5x₃³x₄ + 2.0x₃⁴ + 0.5x₂x₃x₄² + x₂x₃²x₄ - 1.5x₂x₃³ + 0.5x₂²x₃x₄ - 0.5x₂³x₃ + 0.5x₁x₃x₄² + x₁x₃²x₄ - 1.5x₁x₃³ - 1.5x₁x₂x₃x₄ + x₁x₂x₃² + 0.5x₁x₂²x₃ + 0.5x₁²x₃x₄ + 0.5x₁²x₂x₃ - 0.5x₁³x₃ ≥ 0
16: -0.5x₄ + x₄² + 0.5x₃x₄ + 0.5x₂x₄ + 0.5x₁x₄ - 1.5x₄³ + x₃x₄² + 0.5x₃²x₄ + x₂x₄² - 1.5x₂x₃x₄ + 0.5x₂²x₄ + x₁x₄² - 1.5x₁x₃x₄ - 1.5x₁x₂x₄ + 0.5x₁²x₄ + 2.0x₄⁴ - 1.5x₃x₄³ - 0.5x₃³x₄ - 1.5x₂x₄³ + x₂x₃x₄² + 0.5x₂x₃²x₄ + 0.5x₂²x₃x₄ - 0.5x₂³x₄ - 1.5x₁x₄³ + x₁x₃x₄² + 0.5x₁x₃²x₄ + x₁x₂x₄² - 1.5x₁x₂x₃x₄ + 0.5x₁x₂²x₄ + 0.5x₁²x₃x₄ + 0.5x₁²x₂x₄ - 0.5x₁³x₄ ≥ 0"
    end
    if optimize
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, prob, 4).objective ≈ 4 atol = 2e-7
        :COPTMoment ∈ solvers && @test poly_optimize(:COPTMoment, prob, 4).objective ≈ 4 atol = 1e-6
    end
end

@testset "Example 6.7" begin
    DynamicPolynomials.@polyvar x[1:3]
    local prob
    for method in tightening_methods
        prob = roundproblem(x[1]^4 * x[2]^2 + x[2]^4 * x[3]^2 + x[3]^4 * x[1]^2 - 3prod(x .^ 2) + x[2]^2,
            nonneg=[x[1] - x[2]*x[3], -x[2] + x[3]^2], tighter=true)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: x₂² + x₂⁴x₃² + x₁²x₃⁴ - 3.0x₁²x₂²x₃² + x₁⁴x₂²
3 equality constraints
1: 4.0x₂x₃ + 2.0x₂⁴x₃ + 4.0x₁²x₃³ - 6.0x₁²x₂²x₃ + 8.0x₂³x₃³ + 2.0x₁x₂x₃⁴ - 6.0x₁x₂³x₃² - 12.0x₁²x₂x₃³ + 4.0x₁³x₂³ + 4.0x₁⁴x₂x₃ + 4.0x₁x₃⁶ - 12.0x₁x₂²x₃⁴ + 8.0x₁³x₂²x₃² = 0
2: 2.0x₁²x₃⁴ - 6.0x₁²x₂²x₃² + 4.0x₁⁴x₂² - 2.0x₁x₂x₃⁵ + 6.0x₁x₂³x₃³ - 4.0x₁³x₂³x₃ = 0
3: 2.0x₂² - 2.0x₂x₃² + 4.0x₂⁴x₃² - 6.0x₁²x₂²x₃² + 2.0x₁⁴x₂² - 4.0x₂³x₃⁴ + 2.0x₁x₂x₃⁵ - 6.0x₁x₂³x₃³ + 6.0x₁²x₂x₃⁴ + 4.0x₁³x₂³x₃ - 2.0x₁⁴x₂x₃² - 2.0x₁x₃⁷ + 6.0x₁x₂²x₃⁵ - 4.0x₁³x₂²x₃³ = 0
4 nonnegative constraints
4: x₁ - x₂x₃ ≥ 0
5: -x₂ + x₃² ≥ 0
6: 2.0x₁x₃⁴ - 6.0x₁x₂²x₃² + 4.0x₁³x₂² ≥ 0
7: -2.0x₂ - 4.0x₂³x₃² + 6.0x₁²x₂x₃² - 2.0x₁⁴x₂ - 2.0x₁x₃⁵ + 6.0x₁x₂²x₃³ - 4.0x₁³x₂²x₃ ≥ 0"
    end
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test(poly_optimize(solver, prob, 5).objective ≈ 0, atol = solver==:SCSMoment ? 2e-4 : 1e-5,
                    skip = solver ∈ (:ClarabelMoment, :HypatiaMoment))
            end
        end
    end
end

@testset "Example 6.8" begin
    # note that this example is wrong in the paper. The optimum is the smallest root of
    # -2214 + 5697x - 144x^2 + 64x^3 ≈ 0.3918305004
    # and the points are
    # x₁ = x₂ = x₃ = smallest root of -9 + 9x + 4x^3 ≈ 0.7850032632
    # x₄ = smallest root of -125 + 75x + 12x^3 ≈ 1.308338772
    DynamicPolynomials.@polyvar x[1:4]
    local prob
    for method in tightening_methods
        prob = roundproblem(sum(x[i]^2 * (x[i] - x[4])^2 + (x[i] - 1)^2 for i in 1:3) +
            2x[1]*x[2]*x[3]*(x[1] + x[2] + x[3] - 2x[4]), nonneg=[x[1] - x[2], x[2] - x[3]], tighter=method)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 4 variables
Objective: 3.0 - 2.0x₃ - 2.0x₂ - 2.0x₁ + x₃² + x₂² + x₁² + x₃²x₄² - 2.0x₃³x₄ + x₃⁴ + x₂²x₄² - 2.0x₂³x₄ + x₂⁴ - 4.0x₁x₂x₃x₄ + 2.0x₁x₂x₃² + 2.0x₁x₂²x₃ + x₁²x₄² + 2.0x₁²x₂x₃ - 2.0x₁³x₄ + x₁⁴
4 equality constraints
1: -6.0 + 2.0x₃ + 2.0x₂ + 2.0x₁ + 2.0x₃x₄² - 6.0x₃²x₄ + 4.0x₃³ + 2.0x₂x₄² - 4.0x₂x₃x₄ + 2.0x₂x₃² - 6.0x₂²x₄ + 2.0x₂²x₃ + 4.0x₂³ + 2.0x₁x₄² - 4.0x₁x₃x₄ + 2.0x₁x₃² - 4.0x₁x₂x₄ + 12.0x₁x₂x₃ + 2.0x₁x₂² - 6.0x₁²x₄ + 2.0x₁²x₃ + 2.0x₁²x₂ + 4.0x₁³ = 0
2: 2.0x₃²x₄ - 2.0x₃³ + 2.0x₂²x₄ - 2.0x₂³ - 4.0x₁x₂x₃ + 2.0x₁²x₄ - 2.0x₁³ = 0
3: 2.0x₂ - 2.0x₁ - 2.0x₁x₂ + 2.0x₁² + 4.0x₂²x₃x₄ - 2.0x₂²x₃² - 2.0x₂³x₃ - 2.0x₁x₂x₄² - 4.0x₁x₂x₃x₄ + 2.0x₁x₂x₃² - 2.0x₁x₂²x₃ + 2.0x₁²x₄² + 6.0x₁²x₂x₄ + 4.0x₁²x₂x₃ - 6.0x₁³x₄ - 4.0x₁³x₂ + 4.0x₁⁴ = 0
4: -2.0x₃ + 2.0x₂ + 2.0x₃² - 2.0x₂x₃ + 2.0x₃²x₄² - 6.0x₃³x₄ + 4.0x₃⁴ - 2.0x₂x₃x₄² + 6.0x₂x₃²x₄ - 4.0x₂x₃³ - 4.0x₁x₂x₃x₄ + 4.0x₁x₂x₃² + 4.0x₁x₂²x₄ - 2.0x₁x₂²x₃ - 2.0x₁x₂³ + 2.0x₁²x₂x₃ - 2.0x₁²x₂² = 0
4 nonnegative constraints
5: -x₂ + x₁ ≥ 0
6: -x₃ + x₂ ≥ 0
7: -2.0 + 2.0x₁ - 4.0x₂x₃x₄ + 2.0x₂x₃² + 2.0x₂²x₃ + 2.0x₁x₄² + 4.0x₁x₂x₃ - 6.0x₁²x₄ + 4.0x₁³ ≥ 0
8: 2.0 - 2.0x₃ - 2.0x₃x₄² + 6.0x₃²x₄ - 4.0x₃³ + 4.0x₁x₂x₄ - 4.0x₁x₂x₃ - 2.0x₁x₂² - 2.0x₁²x₂ ≥ 0"
    end
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob, 3).objective ≈ 0.3918305004 atol = solver==:SCSMoment ? 2e-3 : 1e-5
            end
        end
    end
end

@testset "Example 6.9" begin
    DynamicPolynomials.@polyvar x[1:4]
    local prob
    for method in tightening_methods
        prob = roundproblem((sum(x) +1)^2 - 4*(x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4] + x[1]), nonneg=[x; 1 .- x],
            tighter=method)
        @test strRep(prob) == "Real-valued polynomial optimization problem in 4 variables
Objective: 1.0 - 2.0x₄ + 2.0x₃ + 2.0x₂ - 2.0x₁ + x₄² - 2.0x₃x₄ + x₃² + 2.0x₂x₄ - 2.0x₂x₃ + x₂² + 2.0x₁x₄ + 2.0x₁x₃ - 2.0x₁x₂ + x₁²
8 equality constraints
01: -2.0x₁ + 2.0x₁x₄ + 2.0x₁x₃ - 2.0x₁x₂ + 4.0x₁² - 2.0x₁²x₄ - 2.0x₁²x₃ + 2.0x₁²x₂ - 2.0x₁³ = 0
02: 2.0x₂ + 2.0x₂x₄ - 2.0x₂x₃ - 2.0x₁x₂ - 2.0x₂²x₄ + 2.0x₂²x₃ - 2.0x₂³ + 2.0x₁x₂² = 0
03: 2.0x₃ - 2.0x₃x₄ - 2.0x₂x₃ + 2.0x₁x₃ + 2.0x₃²x₄ - 2.0x₃³ + 2.0x₂x₃² - 2.0x₁x₃² = 0
04: -2.0x₄ + 4.0x₄² - 2.0x₃x₄ + 2.0x₂x₄ + 2.0x₁x₄ - 2.0x₄³ + 2.0x₃x₄² - 2.0x₂x₄² - 2.0x₁x₄² = 0
05: 2.0x₁ - 2.0x₁x₄ - 2.0x₁x₃ + 2.0x₁x₂ - 4.0x₁² + 2.0x₁²x₄ + 2.0x₁²x₃ - 2.0x₁²x₂ + 2.0x₁³ = 0
06: -2.0x₂ - 2.0x₂x₄ + 2.0x₂x₃ + 2.0x₁x₂ + 2.0x₂²x₄ - 2.0x₂²x₃ + 2.0x₂³ - 2.0x₁x₂² = 0
07: -2.0x₃ + 2.0x₃x₄ + 2.0x₂x₃ - 2.0x₁x₃ - 2.0x₃²x₄ + 2.0x₃³ - 2.0x₂x₃² + 2.0x₁x₃² = 0
08: 2.0x₄ - 4.0x₄² + 2.0x₃x₄ - 2.0x₂x₄ - 2.0x₁x₄ + 2.0x₄³ - 2.0x₃x₄² + 2.0x₂x₄² + 2.0x₁x₄² = 0
16 nonnegative constraints
09: x₁ ≥ 0
10: x₂ ≥ 0
11: x₃ ≥ 0
12: x₄ ≥ 0
13: 1.0 - x₁ ≥ 0
14: 1.0 - x₂ ≥ 0
15: 1.0 - x₃ ≥ 0
16: 1.0 - x₄ ≥ 0
17: -2.0 + 2.0x₄ + 2.0x₃ - 2.0x₂ + 4.0x₁ - 2.0x₁x₄ - 2.0x₁x₃ + 2.0x₁x₂ - 2.0x₁² ≥ 0
18: 2.0 + 2.0x₄ - 2.0x₃ - 2.0x₁ - 2.0x₂x₄ + 2.0x₂x₃ - 2.0x₂² + 2.0x₁x₂ ≥ 0
19: 2.0 - 2.0x₄ - 2.0x₂ + 2.0x₁ + 2.0x₃x₄ - 2.0x₃² + 2.0x₂x₃ - 2.0x₁x₃ ≥ 0
20: -2.0 + 4.0x₄ - 2.0x₃ + 2.0x₂ + 2.0x₁ - 2.0x₄² + 2.0x₃x₄ - 2.0x₂x₄ - 2.0x₁x₄ ≥ 0
21: 2.0x₁ - 2.0x₁x₄ - 2.0x₁x₃ + 2.0x₁x₂ - 2.0x₁² ≥ 0
22: -2.0x₂ - 2.0x₂x₄ + 2.0x₂x₃ - 2.0x₂² + 2.0x₁x₂ ≥ 0
23: -2.0x₃ + 2.0x₃x₄ - 2.0x₃² + 2.0x₂x₃ - 2.0x₁x₃ ≥ 0
24: 2.0x₄ - 2.0x₄² + 2.0x₃x₄ - 2.0x₂x₄ - 2.0x₁x₄ ≥ 0"
    end
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = solver==:SCSMoment ? 1e-3 : 6e-5
            end
        end
    end
end
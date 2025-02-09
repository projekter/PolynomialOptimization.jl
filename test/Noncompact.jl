include("./shared.jl")

@testset "POP 1 (Motzkin)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(x[1]^2 * x[2]^2 * (x[1]^2 + x[2]^2 - 1), noncompact=(1e-5, 2))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 + 5.0e-5x₂² + 5.0e-5x₁² + 0.0001x₂⁴ - 0.9998x₁²x₂² + 0.0001x₁⁴ + 0.0001x₂⁶ - 0.9997x₁²x₂⁴ - 0.9997x₁⁴x₂² + 0.0001x₁⁶ + 5.0e-5x₂⁸ + 1.0002x₁²x₂⁶ + 2.0003x₁⁴x₂⁴ + 1.0002x₁⁶x₂² + 5.0e-5x₁⁸ + 1.0e-5x₂¹⁰ + 1.00005x₁²x₂⁸ + 3.0001x₁⁴x₂⁶ + 3.0001x₁⁶x₂⁴ + 1.00005x₁⁸x₂² + 1.0e-5x₁¹⁰
Objective was scaled by the prefactor 1.0 + 2.0x₂² + 2.0x₁² + x₂⁴ + 2.0x₁²x₂² + x₁⁴"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ -0.0369 atol = 2e-3
            end
        end
    end
end

@testset "POP 2 (Robinson)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(x[1]^6 + x[2]^6 - x[1]^4*x[2]^2 - x[1]^2*x[2]^4 - x[1]^4 - x[2]^4 - x[1]^2 - x[2]^2 + 3x[1]^2*x[2]^2,
        noncompact=(1e-5, 2))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 - 0.99995x₂² - 0.99995x₁² - 2.9999x₂⁴ - 0.9998x₁²x₂² - 2.9999x₁⁴ - 1.9999x₂⁶ + 0.0003x₁²x₂⁴ + 0.0003x₁⁴x₂² - 1.9999x₁⁶ + 1.00005x₂⁸ + 1.0002x₁²x₂⁶ + 0.0003x₁⁴x₂⁴ + 1.0002x₁⁶x₂² + 1.00005x₁⁸ + 1.00001x₂¹⁰ + 1.00005x₁²x₂⁸ - 1.9999x₁⁴x₂⁶ - 1.9999x₁⁶x₂⁴ + 1.00005x₁⁸x₂² + 1.00001x₁¹⁰
Objective was scaled by the prefactor 1.0 + 2.0x₂² + 2.0x₁² + x₂⁴ + 2.0x₁²x₂² + x₁⁴"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ -0.9999 atol = 1e-3
            end
        end
    end
end

@testset "POP 3 (Choi-Lam)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(x[1]^4*x[2]^2 + x[2]^4 + x[1]^2 - 3x[1]^2*x[2]^2, noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 + 4.0e-5x₂² + 1.00004x₁² + 1.00006x₂⁴ - 1.99988x₁²x₂² + 1.00006x₁⁴ + 1.00004x₂⁶ - 1.99988x₁²x₂⁴ - 1.99988x₁⁴x₂² + 4.0e-5x₁⁶ + 1.0e-5x₂⁸ + 4.0e-5x₁²x₂⁶ + 1.00006x₁⁴x₂⁴ + 1.00004x₁⁶x₂² + 1.0e-5x₁⁸
Objective was scaled by the prefactor 1.0 + x₂² + x₁²"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = 3e-4
            end
        end
    end
end

@testset "POP 4 (Lax-Lax)" begin
    DynamicPolynomials.@polyvar x[1:3]
    Random.seed!(1234)
    prob = poly_problem(prod(x) - x[1]*(x[2] - x[1])*(x[3] - x[1])*(1 - x[1]) - x[2]*(x[1] - x[2])*(x[3] - x[2])*(1 - x[2]) -
        x[3]*(x[1] - x[3])*(x[2] - x[3])*(1 - x[3]) - (x[1] -1)*(x[2] -1)*(x[3] -1), noncompact=(1e-5, 2))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: 1.00001 - x₃ - x₂ - x₁ + 2.00004x₃² + x₂x₃ + 2.00004x₂² + x₁x₃ + x₁x₂ + 2.00004x₁² - 3.0x₃³ - x₂x₃² - x₂²x₃ - 3.0x₂³ - x₁x₃² - 3.0x₁x₂x₃ - x₁x₂² - x₁²x₃ - x₁²x₂ - 3.0x₁³ + 2.00006x₃⁴ + x₂x₃³ + 2.00012x₂²x₃² + x₂³x₃ + 2.00006x₂⁴ + x₁x₃³ + 3.0x₁x₂x₃² + 3.0x₁x₂²x₃ + x₁x₂³ + 2.00012x₁²x₃² + 3.0x₁²x₂x₃ + 2.00012x₁²x₂² + x₁³x₃ + x₁³x₂ + 2.00006x₁⁴ - 3.0x₃⁵ + x₂x₃⁴ - 2.0x₂²x₃³ - 2.0x₂³x₃² + x₂⁴x₃ - 3.0x₂⁵ + x₁x₃⁴ - 6.0x₁x₂x₃³ + 2.0x₁x₂²x₃² - 6.0x₁x₂³x₃ + x₁x₂⁴ - 2.0x₁²x₃³ + 2.0x₁²x₂x₃² + 2.0x₁²x₂²x₃ - 2.0x₁²x₂³ - 2.0x₁³x₃² - 6.0x₁³x₂x₃ - 2.0x₁³x₂² + x₁⁴x₃ + x₁⁴x₂ - 3.0x₁⁵ + 2.00004x₃⁶ - x₂x₃⁵ + 2.00012x₂²x₃⁴ - 2.0x₂³x₃³ + 2.00012x₂⁴x₃² - x₂⁵x₃ + 2.00004x₂⁶ - x₁x₃⁵ + 3.0x₁x₂x₃⁴ + 2.0x₁x₂²x₃³ + 2.0x₁x₂³x₃² + 3.0x₁x₂⁴x₃ - x₁x₂⁵ + 2.00012x₁²x₃⁴ + 2.0x₁²x₂x₃³ + 0.00024x₁²x₂²x₃² + 2.0x₁²x₂³x₃ + 2.00012x₁²x₂⁴ - 2.0x₁³x₃³ + 2.0x₁³x₂x₃² + 2.0x₁³x₂²x₃ - 2.0x₁³x₂³ + 2.00012x₁⁴x₃² + 3.0x₁⁴x₂x₃ + 2.00012x₁⁴x₂² - x₁⁵x₃ - x₁⁵x₂ + 2.00004x₁⁶ - x₃⁷ + x₂x₃⁶ - x₂²x₃⁵ + x₂³x₃⁴ + x₂⁴x₃³ - x₂⁵x₃² + x₂⁶x₃ - x₂⁷ + x₁x₃⁶ - 3.0x₁x₂x₃⁵ + 3.0x₁x₂²x₃⁴ - 6.0x₁x₂³x₃³ + 3.0x₁x₂⁴x₃² - 3.0x₁x₂⁵x₃ + x₁x₂⁶ - x₁²x₃⁵ + 3.0x₁²x₂x₃⁴ + 2.0x₁²x₂²x₃³ + 2.0x₁²x₂³x₃² + 3.0x₁²x₂⁴x₃ - x₁²x₂⁵ + x₁³x₃⁴ - 6.0x₁³x₂x₃³ + 2.0x₁³x₂²x₃² - 6.0x₁³x₂³x₃ + x₁³x₂⁴ + x₁⁴x₃³ + 3.0x₁⁴x₂x₃² + 3.0x₁⁴x₂²x₃ + x₁⁴x₂³ - x₁⁵x₃² - 3.0x₁⁵x₂x₃ - x₁⁵x₂² + x₁⁶x₃ + x₁⁶x₂ - x₁⁷ + 1.00001x₃⁸ - x₂x₃⁷ + 2.00004x₂²x₃⁶ - 3.0x₂³x₃⁵ + 2.00006x₂⁴x₃⁴ - 3.0x₂⁵x₃³ + 2.00004x₂⁶x₃² - x₂⁷x₃ + 1.00001x₂⁸ - x₁x₃⁷ + x₁x₂x₃⁶ - x₁x₂²x₃⁵ + x₁x₂³x₃⁴ + x₁x₂⁴x₃³ - x₁x₂⁵x₃² + x₁x₂⁶x₃ - x₁x₂⁷ + 2.00004x₁²x₃⁶ - x₁²x₂x₃⁵ + 2.00012x₁²x₂²x₃⁴ - 2.0x₁²x₂³x₃³ + 2.00012x₁²x₂⁴x₃² - x₁²x₂⁵x₃ + 2.00004x₁²x₂⁶ - 3.0x₁³x₃⁵ + x₁³x₂x₃⁴ - 2.0x₁³x₂²x₃³ - 2.0x₁³x₂³x₃² + x₁³x₂⁴x₃ - 3.0x₁³x₂⁵ + 2.00006x₁⁴x₃⁴ + x₁⁴x₂x₃³ + 2.00012x₁⁴x₂²x₃² + x₁⁴x₂³x₃ + 2.00006x₁⁴x₂⁴ - 3.0x₁⁵x₃³ - x₁⁵x₂x₃² - x₁⁵x₂²x₃ - 3.0x₁⁵x₂³ + 2.00004x₁⁶x₃² + x₁⁶x₂x₃ + 2.00004x₁⁶x₂² - x₁⁷x₃ - x₁⁷x₂ + 1.00001x₁⁸
Objective was scaled by the prefactor 1.0 + 2.0x₃² + 2.0x₂² + 2.0x₁² + x₃⁴ + 2.0x₂²x₃² + x₂⁴ + 2.0x₁²x₃² + 2.0x₁²x₂² + x₁⁴"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = 2e-3
            end
        end
    end
end

@testset "POP 5 (Delzell)" begin
    DynamicPolynomials.@polyvar x[1:3]
    Random.seed!(1234)
    prob = poly_problem(x[1]^4*x[2]^2 + x[2]^4*x[3]^2 + x[1]^2*x[3]^4 - 3x[1]^2*x[2]^2*x[3]^2 + x[3]^8, noncompact=(1e-5, 2))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: 1.0e-5 + 6.0e-5x₃² + 6.0e-5x₂² + 6.0e-5x₁² + 0.00015x₃⁴ + 0.0003x₂²x₃² + 0.00015x₂⁴ + 0.0003x₁²x₃² + 0.0003x₁²x₂² + 0.00015x₁⁴ + 0.0002x₃⁶ + 0.0006x₂²x₃⁴ + 1.0006x₂⁴x₃² + 0.0002x₂⁶ + 1.0006x₁²x₃⁴ - 2.9988x₁²x₂²x₃² + 0.0006x₁²x₂⁴ + 0.0006x₁⁴x₃² + 1.0006x₁⁴x₂² + 0.0002x₁⁶ + 1.00015x₃⁸ + 0.0006x₂²x₃⁶ + 2.0009x₂⁴x₃⁴ + 2.0006x₂⁶x₃² + 0.00015x₂⁸ + 2.0006x₁²x₃⁶ - 3.9982x₁²x₂²x₃⁴ - 3.9982x₁²x₂⁴x₃² + 0.0006x₁²x₂⁶ + 2.0009x₁⁴x₃⁴ - 3.9982x₁⁴x₂²x₃² + 2.0009x₁⁴x₂⁴ + 0.0006x₁⁶x₃² + 2.0006x₁⁶x₂² + 0.00015x₁⁸ + 2.00006x₃¹⁰ + 2.0003x₂²x₃⁸ + 1.0006x₂⁴x₃⁶ + 2.0006x₂⁶x₃⁴ + 1.0003x₂⁸x₃² + 6.0e-5x₂¹⁰ + 3.0003x₁²x₃⁸ - 0.9988x₁²x₂²x₃⁶ - 2.9982x₁²x₂⁴x₃⁴ - 0.9988x₁²x₂⁶x₃² + 0.0003x₁²x₂⁸ + 2.0006x₁⁴x₃⁶ - 2.9982x₁⁴x₂²x₃⁴ - 2.9982x₁⁴x₂⁴x₃² + 1.0006x₁⁴x₂⁶ + 1.0006x₁⁶x₃⁴ - 0.9988x₁⁶x₂²x₃² + 2.0006x₁⁶x₂⁴ + 0.0003x₁⁸x₃² + 1.0003x₁⁸x₂² + 6.0e-5x₁¹⁰ + 1.00001x₃¹² + 2.00006x₂²x₃¹⁰ + 1.00015x₂⁴x₃⁸ + 0.0002x₂⁶x₃⁶ + 0.00015x₂⁸x₃⁴ + 6.0e-5x₂¹⁰x₃² + 1.0e-5x₂¹² + 2.00006x₁²x₃¹⁰ + 2.0003x₁²x₂²x₃⁸ + 0.0006x₁²x₂⁴x₃⁶ + 0.0006x₁²x₂⁶x₃⁴ + 0.0003x₁²x₂⁸x₃² + 6.0e-5x₁²x₂¹⁰ + 1.00015x₁⁴x₃⁸ + 0.0006x₁⁴x₂²x₃⁶ + 0.0009x₁⁴x₂⁴x₃⁴ + 0.0006x₁⁴x₂⁶x₃² + 0.00015x₁⁴x₂⁸ + 0.0002x₁⁶x₃⁶ + 0.0006x₁⁶x₂²x₃⁴ + 0.0006x₁⁶x₂⁴x₃² + 0.0002x₁⁶x₂⁶ + 0.00015x₁⁸x₃⁴ + 0.0003x₁⁸x₂²x₃² + 0.00015x₁⁸x₂⁴ + 6.0e-5x₁¹⁰x₃² + 6.0e-5x₁¹⁰x₂² + 1.0e-5x₁¹²
Objective was scaled by the prefactor 1.0 + 2.0x₃² + 2.0x₂² + 2.0x₁² + x₃⁴ + 2.0x₂²x₃² + x₂⁴ + 2.0x₁²x₃² + 2.0x₁²x₂² + x₁⁴"
    if optimize
       for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = 1e-3
            end
       end
    end
end

@testset "POP 6 (Modified Motzkin)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem((x[1]^2 + x[2]^2 - 3)*x[1]^2*x[2]^2, nonneg=[x[1]^2+x[2]^2-4], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 + 5.0e-5x₂² + 5.0e-5x₁² + 0.0001x₂⁴ - 2.9998x₁²x₂² + 0.0001x₁⁴ + 0.0001x₂⁶ - 1.9997x₁²x₂⁴ - 1.9997x₁⁴x₂² + 0.0001x₁⁶ + 5.0e-5x₂⁸ + 1.0002x₁²x₂⁶ + 2.0003x₁⁴x₂⁴ + 1.0002x₁⁶x₂² + 5.0e-5x₁⁸ + 1.0e-5x₂¹⁰ + 5.0e-5x₁²x₂⁸ + 0.0001x₁⁴x₂⁶ + 0.0001x₁⁶x₂⁴ + 5.0e-5x₁⁸x₂² + 1.0e-5x₁¹⁰
Objective was scaled by the prefactor 1.0 + x₂² + x₁²
1 nonnegative constraint
1: -4.0 + x₂² + x₁² ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0.0062 atol = 2e-2
            end
        end
    end
end

@testset "POP 7 (Example 4.3)" begin
    DynamicPolynomials.@polyvar x[1:3]
    Random.seed!(1234)
    prob = poly_problem(x[1]^4 + x[2]^4 + x[3]^4 - 4.0x[1]*x[3]^3, nonneg=[1-x[1]^4+.5x[2]^4-x[3]^4], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: 1.0e-5 + 4.0e-5x₃² + 4.0e-5x₂² + 4.0e-5x₁² + 1.00006x₃⁴ + 0.00012x₂²x₃² + 1.00006x₂⁴ - 4.0x₁x₃³ + 0.00012x₁²x₃² + 0.00012x₁²x₂² + 1.00006x₁⁴ + 1.00004x₃⁶ + 1.00012x₂²x₃⁴ + 1.00012x₂⁴x₃² + 1.00004x₂⁶ - 4.0x₁x₃⁵ - 4.0x₁x₂²x₃³ + 1.00012x₁²x₃⁴ + 0.00024x₁²x₂²x₃² + 1.00012x₁²x₂⁴ - 4.0x₁³x₃³ + 1.00012x₁⁴x₃² + 1.00012x₁⁴x₂² + 1.00004x₁⁶ + 1.0e-5x₃⁸ + 4.0e-5x₂²x₃⁶ + 6.0e-5x₂⁴x₃⁴ + 4.0e-5x₂⁶x₃² + 1.0e-5x₂⁸ + 4.0e-5x₁²x₃⁶ + 0.00012x₁²x₂²x₃⁴ + 0.00012x₁²x₂⁴x₃² + 4.0e-5x₁²x₂⁶ + 6.0e-5x₁⁴x₃⁴ + 0.00012x₁⁴x₂²x₃² + 6.0e-5x₁⁴x₂⁴ + 4.0e-5x₁⁶x₃² + 4.0e-5x₁⁶x₂² + 1.0e-5x₁⁸
Objective was scaled by the prefactor 1.0 + x₃² + x₂² + x₁²
1 nonnegative constraint
1: 1.0 - x₃⁴ + 0.5x₂⁴ - x₁⁴ ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ -1.27937458 atol = solver==:SCSMoment ? 1e-2 : 1e-7
            end
        end
    end
end

@testset "POP 8 (Example 3.1)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(x[1]^2 + 1.0, nonneg=[1-x[2]^2, x[2]^2-0.25], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.00001 + 1.00003x₂² + 2.00003x₁² + 3.0e-5x₂⁴ + 1.00006x₁²x₂² + 1.00003x₁⁴ + 1.0e-5x₂⁶ + 3.0e-5x₁²x₂⁴ + 3.0e-5x₁⁴x₂² + 1.0e-5x₁⁶
Objective was scaled by the prefactor 1.0 + x₂² + x₁²
2 nonnegative constraints
1: 1.0 - x₂² ≥ 0
2: -0.25 + x₂² ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 1 atol = 1e-4
            end
        end
    end
end

@testset "POP 9 (Example 4.5)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(x[1]^2 + x[2]^2, nonneg=[x[1]^2-x[1]*x[2]-1, x[1]^2+x[1]*x[2]-1, x[2]^2-1], noncompact=(1e-5, 2))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 + 1.00004x₂² + 1.00004x₁² + 2.00006x₂⁴ + 4.00012x₁²x₂² + 2.00006x₁⁴ + 1.00004x₂⁶ + 3.00012x₁²x₂⁴ + 3.00012x₁⁴x₂² + 1.00004x₁⁶ + 1.0e-5x₂⁸ + 4.0e-5x₁²x₂⁶ + 6.0e-5x₁⁴x₂⁴ + 4.0e-5x₁⁶x₂² + 1.0e-5x₁⁸
Objective was scaled by the prefactor 1.0 + 2.0x₂² + 2.0x₁² + x₂⁴ + 2.0x₁²x₂² + x₁⁴
3 nonnegative constraints
1: -1.0 - x₁x₂ + x₁² ≥ 0
2: -1.0 + x₁x₂ + x₁² ≥ 0
3: -1.0 + x₂² ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 3.6182472 atol = 2e-4
            end
        end
    end
end

@testset "POP 10 (Example 4.4)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(-4x[1]^2/3+2x[2]^2/3-2x[1]*x[2], nonneg=[x[2]^2-x[1]^2, -x[1]*x[2]], noncompact=(1e-5, 5))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 + 0.66674x₂² - 2.0x₁x₂ - 1.33326x₁² + 3.33354x₂⁴ - 10.0x₁x₂³ - 3.33291x₁²x₂² - 10.0x₁³x₂ - 6.66646x₁⁴ + 6.66702x₂⁶ - 20.0x₁x₂⁵ + 0.00105x₁²x₂⁴ - 40.0x₁³x₂³ - 19.99895x₁⁴x₂² - 20.0x₁⁵x₂ - 13.33298x₁⁶ + 6.66702x₂⁸ - 20.0x₁x₂⁷ + 6.66807x₁²x₂⁶ - 60.0x₁³x₂⁵ - 19.9979x₁⁴x₂⁴ - 60.0x₁⁵x₂³ - 33.33193x₁⁶x₂² - 20.0x₁⁷x₂ - 13.33298x₁⁸ + 3.33354x₂¹⁰ - 10.0x₁x₂⁹ + 6.66772x₁²x₂⁸ - 40.0x₁³x₂⁷ - 6.66457x₁⁴x₂⁶ - 60.0x₁⁵x₂⁵ - 26.66457x₁⁶x₂⁴ - 40.0x₁⁷x₂³ - 23.33228x₁⁸x₂² - 10.0x₁⁹x₂ - 6.66646x₁¹⁰ + 0.66674x₂¹² - 2.0x₁x₂¹¹ + 2.00042x₁²x₂¹⁰ - 10.0x₁³x₂⁹ + 0.00105x₁⁴x₂⁸ - 20.0x₁⁵x₂⁷ - 6.66527x₁⁶x₂⁶ - 20.0x₁⁷x₂⁵ - 9.99895x₁⁸x₂⁴ - 10.0x₁⁹x₂³ - 5.99958x₁¹⁰x₂² - 2.0x₁¹¹x₂ - 1.33326x₁¹² + 1.0e-5x₂¹⁴ + 7.0e-5x₁²x₂¹² + 0.00021x₁⁴x₂¹⁰ + 0.00035x₁⁶x₂⁸ + 0.00035x₁⁸x₂⁶ + 0.00021x₁¹⁰x₂⁴ + 7.0e-5x₁¹²x₂² + 1.0e-5x₁¹⁴
Objective was scaled by the prefactor 1.0 + 5.0x₂² + 5.0x₁² + 10.0x₂⁴ + 20.0x₁²x₂² + 10.0x₁⁴ + 10.0x₂⁶ + 30.0x₁²x₂⁴ + 30.0x₁⁴x₂² + 10.0x₁⁶ + 5.0x₂⁸ + 20.0x₁²x₂⁶ + 30.0x₁⁴x₂⁴ + 20.0x₁⁶x₂² + 5.0x₁⁸ + x₂¹⁰ + 5.0x₁²x₂⁸ + 10.0x₁⁴x₂⁶ + 10.0x₁⁶x₂⁴ + 5.0x₁⁸x₂² + x₁¹⁰
2 nonnegative constraints
1: x₂² - x₁² ≥ 0
2: -x₁x₂ ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = 2e-2
            end
        end
    end
end

@testset "POP 11 (§5.2)" begin
    # Note that the value 1.0072 reported in the document is wrong. Using their original code, we also get -112.
    DynamicPolynomials.@polyvar x[1:8]
    Random.seed!(1234)
    prob = poly_problem(1 + sum((x[j] - x[j-1]^2)^2 + (1-x[j]^2) for j in 2:8), nonneg=x, noncompact=(1e-5, 0))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 8 variables
Objective: 8.00001 + 3.0e-5x₈² + 3.0e-5x₇² + 3.0e-5x₆² + 3.0e-5x₅² + 3.0e-5x₄² + 3.0e-5x₃² + 3.0e-5x₂² + 3.0e-5x₁² - 2.0x₇²x₈ - 2.0x₆²x₇ - 2.0x₅²x₆ - 2.0x₄²x₅ - 2.0x₃²x₄ - 2.0x₂²x₃ - 2.0x₁²x₂ + 3.0e-5x₈⁴ + 6.0e-5x₇²x₈² + 1.00003x₇⁴ + 6.0e-5x₆²x₈² + 6.0e-5x₆²x₇² + 1.00003x₆⁴ + 6.0e-5x₅²x₈² + 6.0e-5x₅²x₇² + 6.0e-5x₅²x₆² + 1.00003x₅⁴ + 6.0e-5x₄²x₈² + 6.0e-5x₄²x₇² + 6.0e-5x₄²x₆² + 6.0e-5x₄²x₅² + 1.00003x₄⁴ + 6.0e-5x₃²x₈² + 6.0e-5x₃²x₇² + 6.0e-5x₃²x₆² + 6.0e-5x₃²x₅² + 6.0e-5x₃²x₄² + 1.00003x₃⁴ + 6.0e-5x₂²x₈² + 6.0e-5x₂²x₇² + 6.0e-5x₂²x₆² + 6.0e-5x₂²x₅² + 6.0e-5x₂²x₄² + 6.0e-5x₂²x₃² + 1.00003x₂⁴ + 6.0e-5x₁²x₈² + 6.0e-5x₁²x₇² + 6.0e-5x₁²x₆² + 6.0e-5x₁²x₅² + 6.0e-5x₁²x₄² + 6.0e-5x₁²x₃² + 6.0e-5x₁²x₂² + 1.00003x₁⁴ + 1.0e-5x₈⁶ + 3.0e-5x₇²x₈⁴ + 3.0e-5x₇⁴x₈² + 1.0e-5x₇⁶ + 3.0e-5x₆²x₈⁴ + 6.0e-5x₆²x₇²x₈² + 3.0e-5x₆²x₇⁴ + 3.0e-5x₆⁴x₈² + 3.0e-5x₆⁴x₇² + 1.0e-5x₆⁶ + 3.0e-5x₅²x₈⁴ + 6.0e-5x₅²x₇²x₈² + 3.0e-5x₅²x₇⁴ + 6.0e-5x₅²x₆²x₈² + 6.0e-5x₅²x₆²x₇² + 3.0e-5x₅²x₆⁴ + 3.0e-5x₅⁴x₈² + 3.0e-5x₅⁴x₇² + 3.0e-5x₅⁴x₆² + 1.0e-5x₅⁶ + 3.0e-5x₄²x₈⁴ + 6.0e-5x₄²x₇²x₈² + 3.0e-5x₄²x₇⁴ + 6.0e-5x₄²x₆²x₈² + 6.0e-5x₄²x₆²x₇² + 3.0e-5x₄²x₆⁴ + 6.0e-5x₄²x₅²x₈² + 6.0e-5x₄²x₅²x₇² + 6.0e-5x₄²x₅²x₆² + 3.0e-5x₄²x₅⁴ + 3.0e-5x₄⁴x₈² + 3.0e-5x₄⁴x₇² + 3.0e-5x₄⁴x₆² + 3.0e-5x₄⁴x₅² + 1.0e-5x₄⁶ + 3.0e-5x₃²x₈⁴ + 6.0e-5x₃²x₇²x₈² + 3.0e-5x₃²x₇⁴ + 6.0e-5x₃²x₆²x₈² + 6.0e-5x₃²x₆²x₇² + 3.0e-5x₃²x₆⁴ + 6.0e-5x₃²x₅²x₈² + 6.0e-5x₃²x₅²x₇² + 6.0e-5x₃²x₅²x₆² + 3.0e-5x₃²x₅⁴ + 6.0e-5x₃²x₄²x₈² + 6.0e-5x₃²x₄²x₇² + 6.0e-5x₃²x₄²x₆² + 6.0e-5x₃²x₄²x₅² + 3.0e-5x₃²x₄⁴ + 3.0e-5x₃⁴x₈² + 3.0e-5x₃⁴x₇² + 3.0e-5x₃⁴x₆² + 3.0e-5x₃⁴x₅² + 3.0e-5x₃⁴x₄² + 1.0e-5x₃⁶ + 3.0e-5x₂²x₈⁴ + 6.0e-5x₂²x₇²x₈² + 3.0e-5x₂²x₇⁴ + 6.0e-5x₂²x₆²x₈² + 6.0e-5x₂²x₆²x₇² + 3.0e-5x₂²x₆⁴ + 6.0e-5x₂²x₅²x₈² + 6.0e-5x₂²x₅²x₇² + 6.0e-5x₂²x₅²x₆² + 3.0e-5x₂²x₅⁴ + 6.0e-5x₂²x₄²x₈² + 6.0e-5x₂²x₄²x₇² + 6.0e-5x₂²x₄²x₆² + 6.0e-5x₂²x₄²x₅² + 3.0e-5x₂²x₄⁴ + 6.0e-5x₂²x₃²x₈² + 6.0e-5x₂²x₃²x₇² + 6.0e-5x₂²x₃²x₆² + 6.0e-5x₂²x₃²x₅² + 6.0e-5x₂²x₃²x₄² + 3.0e-5x₂²x₃⁴ + 3.0e-5x₂⁴x₈² + 3.0e-5x₂⁴x₇² + 3.0e-5x₂⁴x₆² + 3.0e-5x₂⁴x₅² + 3.0e-5x₂⁴x₄² + 3.0e-5x₂⁴x₃² + 1.0e-5x₂⁶ + 3.0e-5x₁²x₈⁴ + 6.0e-5x₁²x₇²x₈² + 3.0e-5x₁²x₇⁴ + 6.0e-5x₁²x₆²x₈² + 6.0e-5x₁²x₆²x₇² + 3.0e-5x₁²x₆⁴ + 6.0e-5x₁²x₅²x₈² + 6.0e-5x₁²x₅²x₇² + 6.0e-5x₁²x₅²x₆² + 3.0e-5x₁²x₅⁴ + 6.0e-5x₁²x₄²x₈² + 6.0e-5x₁²x₄²x₇² + 6.0e-5x₁²x₄²x₆² + 6.0e-5x₁²x₄²x₅² + 3.0e-5x₁²x₄⁴ + 6.0e-5x₁²x₃²x₈² + 6.0e-5x₁²x₃²x₇² + 6.0e-5x₁²x₃²x₆² + 6.0e-5x₁²x₃²x₅² + 6.0e-5x₁²x₃²x₄² + 3.0e-5x₁²x₃⁴ + 6.0e-5x₁²x₂²x₈² + 6.0e-5x₁²x₂²x₇² + 6.0e-5x₁²x₂²x₆² + 6.0e-5x₁²x₂²x₅² + 6.0e-5x₁²x₂²x₄² + 6.0e-5x₁²x₂²x₃² + 3.0e-5x₁²x₂⁴ + 3.0e-5x₁⁴x₈² + 3.0e-5x₁⁴x₇² + 3.0e-5x₁⁴x₆² + 3.0e-5x₁⁴x₅² + 3.0e-5x₁⁴x₄² + 3.0e-5x₁⁴x₃² + 3.0e-5x₁⁴x₂² + 1.0e-5x₁⁶
8 nonnegative constraints
1: x₁ ≥ 0
2: x₂ ≥ 0
3: x₃ ≥ 0
4: x₄ ≥ 0
5: x₅ ≥ 0
6: x₆ ≥ 0
7: x₇ ≥ 0
8: x₈ ≥ 0"
    # this problem reformulation seems to be close to ill-posed
    if optimize
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, prob).objective ≈ -112.014 atol = 2e-3
        # COPT: imprecise, Clarabel: no chance with default solver/extremely long with Pardiso, MosekMoment: unknown,
        # HypatiaMoment (dense, sparse OOM): 6min, SCS: very slow
    end
end

@testset "POP 12 (§5.3)" begin
    DynamicPolynomials.@polyvar x[1:8]
    Random.seed!(1234)
    prob = poly_problem(1 + sum((x[2l] - x[2l-1]^2)^2 +
                                (1 - x[2l-1])^2 +
                                90(x[2l+2]^2 - x[2l+1])^2 +
                                (x[2l+1] -1)^2 +
                                10(x[2l] + x[2l+2] -2)^2 +
                                .1(x[2l] - x[2l+2])^2 for l in 1:3), nonneg=x, noncompact=(1e-5, 0))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 8 variables
Objective: 127.00001 - 40.0x₈ - 2.0x₇ - 80.0x₆ - 4.0x₅ - 80.0x₄ - 4.0x₃ - 40.0x₂ - 2.0x₁ + 10.10003x₈² + 91.00003x₇² + 19.8x₆x₈ + 21.20003x₆² + 92.00003x₅² + 19.8x₄x₆ + 21.20003x₄² + 92.00003x₃² + 19.8x₂x₄ + 11.10003x₂² + 1.00003x₁² - 180.0x₇x₈² - 180.0x₅x₆² - 2.0x₅²x₆ - 180.0x₃x₄² - 2.0x₃²x₄ - 2.0x₁²x₂ + 90.00003x₈⁴ + 6.0e-5x₇²x₈² + 3.0e-5x₇⁴ + 6.0e-5x₆²x₈² + 6.0e-5x₆²x₇² + 90.00003x₆⁴ + 6.0e-5x₅²x₈² + 6.0e-5x₅²x₇² + 6.0e-5x₅²x₆² + 1.00003x₅⁴ + 6.0e-5x₄²x₈² + 6.0e-5x₄²x₇² + 6.0e-5x₄²x₆² + 6.0e-5x₄²x₅² + 90.00003x₄⁴ + 6.0e-5x₃²x₈² + 6.0e-5x₃²x₇² + 6.0e-5x₃²x₆² + 6.0e-5x₃²x₅² + 6.0e-5x₃²x₄² + 1.00003x₃⁴ + 6.0e-5x₂²x₈² + 6.0e-5x₂²x₇² + 6.0e-5x₂²x₆² + 6.0e-5x₂²x₅² + 6.0e-5x₂²x₄² + 6.0e-5x₂²x₃² + 3.0e-5x₂⁴ + 6.0e-5x₁²x₈² + 6.0e-5x₁²x₇² + 6.0e-5x₁²x₆² + 6.0e-5x₁²x₅² + 6.0e-5x₁²x₄² + 6.0e-5x₁²x₃² + 6.0e-5x₁²x₂² + 1.00003x₁⁴ + 1.0e-5x₈⁶ + 3.0e-5x₇²x₈⁴ + 3.0e-5x₇⁴x₈² + 1.0e-5x₇⁶ + 3.0e-5x₆²x₈⁴ + 6.0e-5x₆²x₇²x₈² + 3.0e-5x₆²x₇⁴ + 3.0e-5x₆⁴x₈² + 3.0e-5x₆⁴x₇² + 1.0e-5x₆⁶ + 3.0e-5x₅²x₈⁴ + 6.0e-5x₅²x₇²x₈² + 3.0e-5x₅²x₇⁴ + 6.0e-5x₅²x₆²x₈² + 6.0e-5x₅²x₆²x₇² + 3.0e-5x₅²x₆⁴ + 3.0e-5x₅⁴x₈² + 3.0e-5x₅⁴x₇² + 3.0e-5x₅⁴x₆² + 1.0e-5x₅⁶ + 3.0e-5x₄²x₈⁴ + 6.0e-5x₄²x₇²x₈² + 3.0e-5x₄²x₇⁴ + 6.0e-5x₄²x₆²x₈² + 6.0e-5x₄²x₆²x₇² + 3.0e-5x₄²x₆⁴ + 6.0e-5x₄²x₅²x₈² + 6.0e-5x₄²x₅²x₇² + 6.0e-5x₄²x₅²x₆² + 3.0e-5x₄²x₅⁴ + 3.0e-5x₄⁴x₈² + 3.0e-5x₄⁴x₇² + 3.0e-5x₄⁴x₆² + 3.0e-5x₄⁴x₅² + 1.0e-5x₄⁶ + 3.0e-5x₃²x₈⁴ + 6.0e-5x₃²x₇²x₈² + 3.0e-5x₃²x₇⁴ + 6.0e-5x₃²x₆²x₈² + 6.0e-5x₃²x₆²x₇² + 3.0e-5x₃²x₆⁴ + 6.0e-5x₃²x₅²x₈² + 6.0e-5x₃²x₅²x₇² + 6.0e-5x₃²x₅²x₆² + 3.0e-5x₃²x₅⁴ + 6.0e-5x₃²x₄²x₈² + 6.0e-5x₃²x₄²x₇² + 6.0e-5x₃²x₄²x₆² + 6.0e-5x₃²x₄²x₅² + 3.0e-5x₃²x₄⁴ + 3.0e-5x₃⁴x₈² + 3.0e-5x₃⁴x₇² + 3.0e-5x₃⁴x₆² + 3.0e-5x₃⁴x₅² + 3.0e-5x₃⁴x₄² + 1.0e-5x₃⁶ + 3.0e-5x₂²x₈⁴ + 6.0e-5x₂²x₇²x₈² + 3.0e-5x₂²x₇⁴ + 6.0e-5x₂²x₆²x₈² + 6.0e-5x₂²x₆²x₇² + 3.0e-5x₂²x₆⁴ + 6.0e-5x₂²x₅²x₈² + 6.0e-5x₂²x₅²x₇² + 6.0e-5x₂²x₅²x₆² + 3.0e-5x₂²x₅⁴ + 6.0e-5x₂²x₄²x₈² + 6.0e-5x₂²x₄²x₇² + 6.0e-5x₂²x₄²x₆² + 6.0e-5x₂²x₄²x₅² + 3.0e-5x₂²x₄⁴ + 6.0e-5x₂²x₃²x₈² + 6.0e-5x₂²x₃²x₇² + 6.0e-5x₂²x₃²x₆² + 6.0e-5x₂²x₃²x₅² + 6.0e-5x₂²x₃²x₄² + 3.0e-5x₂²x₃⁴ + 3.0e-5x₂⁴x₈² + 3.0e-5x₂⁴x₇² + 3.0e-5x₂⁴x₆² + 3.0e-5x₂⁴x₅² + 3.0e-5x₂⁴x₄² + 3.0e-5x₂⁴x₃² + 1.0e-5x₂⁶ + 3.0e-5x₁²x₈⁴ + 6.0e-5x₁²x₇²x₈² + 3.0e-5x₁²x₇⁴ + 6.0e-5x₁²x₆²x₈² + 6.0e-5x₁²x₆²x₇² + 3.0e-5x₁²x₆⁴ + 6.0e-5x₁²x₅²x₈² + 6.0e-5x₁²x₅²x₇² + 6.0e-5x₁²x₅²x₆² + 3.0e-5x₁²x₅⁴ + 6.0e-5x₁²x₄²x₈² + 6.0e-5x₁²x₄²x₇² + 6.0e-5x₁²x₄²x₆² + 6.0e-5x₁²x₄²x₅² + 3.0e-5x₁²x₄⁴ + 6.0e-5x₁²x₃²x₈² + 6.0e-5x₁²x₃²x₇² + 6.0e-5x₁²x₃²x₆² + 6.0e-5x₁²x₃²x₅² + 6.0e-5x₁²x₃²x₄² + 3.0e-5x₁²x₃⁴ + 6.0e-5x₁²x₂²x₈² + 6.0e-5x₁²x₂²x₇² + 6.0e-5x₁²x₂²x₆² + 6.0e-5x₁²x₂²x₅² + 6.0e-5x₁²x₂²x₄² + 6.0e-5x₁²x₂²x₃² + 3.0e-5x₁²x₂⁴ + 3.0e-5x₁⁴x₈² + 3.0e-5x₁⁴x₇² + 3.0e-5x₁⁴x₆² + 3.0e-5x₁⁴x₅² + 3.0e-5x₁⁴x₄² + 3.0e-5x₁⁴x₃² + 3.0e-5x₁⁴x₂² + 1.0e-5x₁⁶
8 nonnegative constraints
1: x₁ ≥ 0
2: x₂ ≥ 0
3: x₃ ≥ 0
4: x₄ ≥ 0
5: x₅ ≥ 0
6: x₆ ≥ 0
7: x₇ ≥ 0
8: x₈ ≥ 0"
    if optimize
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, prob).objective ≈ 1.00728 atol = 1e-5
        :COPTMoment ∈ solvers && @test poly_optimize(:COPTMoment, prob).objective ≈ 1.00728 atol = 1e-5
        # MosekMoment: 5min, HypatiaMoment (dense): 4min, Clarabel (Pardiso): 26min
    end
end

@testset "POP 13 (Example A.2)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem((x[1]^2 + x[2]^2 -2)*(x[1]^2 + x[2]^2), zero=[(x[1]^2+x[2]^2-1)*(x[1]-3)], noncompact=(1e-5, 2))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 - 1.99995x₂² - 1.99995x₁² - 2.9999x₂⁴ - 5.9998x₁²x₂² - 2.9999x₁⁴ + 0.0001x₂⁶ + 0.0003x₁²x₂⁴ + 0.0003x₁⁴x₂² + 0.0001x₁⁶ + 1.00005x₂⁸ + 4.0002x₁²x₂⁶ + 6.0003x₁⁴x₂⁴ + 4.0002x₁⁶x₂² + 1.00005x₁⁸ + 1.0e-5x₂¹⁰ + 5.0e-5x₁²x₂⁸ + 0.0001x₁⁴x₂⁶ + 0.0001x₁⁶x₂⁴ + 5.0e-5x₁⁸x₂² + 1.0e-5x₁¹⁰
Objective was scaled by the prefactor 1.0 + 2.0x₂² + 2.0x₁² + x₂⁴ + 2.0x₁²x₂² + x₁⁴
1 equality constraint
1: 3.0 - x₁ - 3.0x₂² - 3.0x₁² + x₁x₂² + x₁³ = 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ -1 atol = 1e-3
            end
        end
    end
end

@testset "POP 14 (Example A.5)" begin
    DynamicPolynomials.@polyvar x[1:3]
    Random.seed!(1234)
    prob = poly_problem(sum(x .^ 6) + 3prod(x .^ 2) - x[1]^2*(x[2]^4 + x[3]^4) - x[2]^2*(x[3]^4 + x[1]^4) -
        x[3]^2*(x[1]^4 + x[2]^4), zero=[sum(x)-1], noncompact=(1e-5, 2))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: 1.0e-5 + 6.0e-5x₃² + 6.0e-5x₂² + 6.0e-5x₁² + 0.00015x₃⁴ + 0.0003x₂²x₃² + 0.00015x₂⁴ + 0.0003x₁²x₃² + 0.0003x₁²x₂² + 0.00015x₁⁴ + 1.0002x₃⁶ - 0.9994x₂²x₃⁴ - 0.9994x₂⁴x₃² + 1.0002x₂⁶ - 0.9994x₁²x₃⁴ + 3.0012x₁²x₂²x₃² - 0.9994x₁²x₂⁴ - 0.9994x₁⁴x₃² - 0.9994x₁⁴x₂² + 1.0002x₁⁶ + 2.00015x₃⁸ + 0.0006x₂²x₃⁶ - 3.9991x₂⁴x₃⁴ + 0.0006x₂⁶x₃² + 2.00015x₂⁸ + 0.0006x₁²x₃⁶ + 2.0018x₁²x₂²x₃⁴ + 2.0018x₁²x₂⁴x₃² + 0.0006x₁²x₂⁶ - 3.9991x₁⁴x₃⁴ + 2.0018x₁⁴x₂²x₃² - 3.9991x₁⁴x₂⁴ + 0.0006x₁⁶x₃² + 0.0006x₁⁶x₂² + 2.00015x₁⁸ + 1.00006x₃¹⁰ + 1.0003x₂²x₃⁸ - 1.9994x₂⁴x₃⁶ - 1.9994x₂⁶x₃⁴ + 1.0003x₂⁸x₃² + 1.00006x₂¹⁰ + 1.0003x₁²x₃⁸ + 1.0012x₁²x₂²x₃⁶ + 0.0018x₁²x₂⁴x₃⁴ + 1.0012x₁²x₂⁶x₃² + 1.0003x₁²x₂⁸ - 1.9994x₁⁴x₃⁶ + 0.0018x₁⁴x₂²x₃⁴ + 0.0018x₁⁴x₂⁴x₃² - 1.9994x₁⁴x₂⁶ - 1.9994x₁⁶x₃⁴ + 1.0012x₁⁶x₂²x₃² - 1.9994x₁⁶x₂⁴ + 1.0003x₁⁸x₃² + 1.0003x₁⁸x₂² + 1.00006x₁¹⁰ + 1.0e-5x₃¹² + 6.0e-5x₂²x₃¹⁰ + 0.00015x₂⁴x₃⁸ + 0.0002x₂⁶x₃⁶ + 0.00015x₂⁸x₃⁴ + 6.0e-5x₂¹⁰x₃² + 1.0e-5x₂¹² + 6.0e-5x₁²x₃¹⁰ + 0.0003x₁²x₂²x₃⁸ + 0.0006x₁²x₂⁴x₃⁶ + 0.0006x₁²x₂⁶x₃⁴ + 0.0003x₁²x₂⁸x₃² + 6.0e-5x₁²x₂¹⁰ + 0.00015x₁⁴x₃⁸ + 0.0006x₁⁴x₂²x₃⁶ + 0.0009x₁⁴x₂⁴x₃⁴ + 0.0006x₁⁴x₂⁶x₃² + 0.00015x₁⁴x₂⁸ + 0.0002x₁⁶x₃⁶ + 0.0006x₁⁶x₂²x₃⁴ + 0.0006x₁⁶x₂⁴x₃² + 0.0002x₁⁶x₂⁶ + 0.00015x₁⁸x₃⁴ + 0.0003x₁⁸x₂²x₃² + 0.00015x₁⁸x₂⁴ + 6.0e-5x₁¹⁰x₃² + 6.0e-5x₁¹⁰x₂² + 1.0e-5x₁¹²
Objective was scaled by the prefactor 1.0 + 2.0x₃² + 2.0x₂² + 2.0x₁² + x₃⁴ + 2.0x₂²x₃² + x₂⁴ + 2.0x₁²x₃² + 2.0x₁²x₂² + x₁⁴
1 equality constraint
1: -1.0 + x₃ + x₂ + x₁ = 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = solver==:SCSMoment ? 1e-3 : 5e-5
            end
        end
    end
end

@testset "POP 15 (Example A.6)" begin
    DynamicPolynomials.@polyvar x[1:4]
    Random.seed!(1234)
    prob = poly_problem(prod(x) - sum(x[i]*prod(x[j] - x[i] for j in 1:4 if j != i) for i in 1:4),
        zero=[x[1], x[2]-x[3], x[3]-x[4]], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 4 variables
Objective: 1.0e-5 + 4.0e-5x₄² + 4.0e-5x₃² + 4.0e-5x₂² + 4.0e-5x₁² + 1.00006x₄⁴ - x₃x₄³ + 0.00012x₃²x₄² - x₃³x₄ + 1.00006x₃⁴ - x₂x₄³ + x₂x₃x₄² + x₂x₃²x₄ - x₂x₃³ + 0.00012x₂²x₄² + x₂²x₃x₄ + 0.00012x₂²x₃² - x₂³x₄ - x₂³x₃ + 1.00006x₂⁴ - x₁x₄³ + x₁x₃x₄² + x₁x₃²x₄ - x₁x₃³ + x₁x₂x₄² - 3.0x₁x₂x₃x₄ + x₁x₂x₃² + x₁x₂²x₄ + x₁x₂²x₃ - x₁x₂³ + 0.00012x₁²x₄² + x₁²x₃x₄ + 0.00012x₁²x₃² + x₁²x₂x₄ + x₁²x₂x₃ + 0.00012x₁²x₂² - x₁³x₄ - x₁³x₃ - x₁³x₂ + 1.00006x₁⁴ + 1.00004x₄⁶ - x₃x₄⁵ + 1.00012x₃²x₄⁴ - 2.0x₃³x₄³ + 1.00012x₃⁴x₄² - x₃⁵x₄ + 1.00004x₃⁶ - x₂x₄⁵ + x₂x₃x₄⁴ + x₂x₃⁴x₄ - x₂x₃⁵ + 1.00012x₂²x₄⁴ + 0.00024x₂²x₃²x₄² + 1.00012x₂²x₃⁴ - 2.0x₂³x₄³ - 2.0x₂³x₃³ + 1.00012x₂⁴x₄² + x₂⁴x₃x₄ + 1.00012x₂⁴x₃² - x₂⁵x₄ - x₂⁵x₃ + 1.00004x₂⁶ - x₁x₄⁵ + x₁x₃x₄⁴ + x₁x₃⁴x₄ - x₁x₃⁵ + x₁x₂x₄⁴ - 3.0x₁x₂x₃x₄³ + 2.0x₁x₂x₃²x₄² - 3.0x₁x₂x₃³x₄ + x₁x₂x₃⁴ + 2.0x₁x₂²x₃x₄² + 2.0x₁x₂²x₃²x₄ - 3.0x₁x₂³x₃x₄ + x₁x₂⁴x₄ + x₁x₂⁴x₃ - x₁x₂⁵ + 1.00012x₁²x₄⁴ + 0.00024x₁²x₃²x₄² + 1.00012x₁²x₃⁴ + 2.0x₁²x₂x₃x₄² + 2.0x₁²x₂x₃²x₄ + 0.00024x₁²x₂²x₄² + 2.0x₁²x₂²x₃x₄ + 0.00024x₁²x₂²x₃² + 1.00012x₁²x₂⁴ - 2.0x₁³x₄³ - 2.0x₁³x₃³ - 3.0x₁³x₂x₃x₄ - 2.0x₁³x₂³ + 1.00012x₁⁴x₄² + x₁⁴x₃x₄ + 1.00012x₁⁴x₃² + x₁⁴x₂x₄ + x₁⁴x₂x₃ + 1.00012x₁⁴x₂² - x₁⁵x₄ - x₁⁵x₃ - x₁⁵x₂ + 1.00004x₁⁶ + 1.0e-5x₄⁸ + 4.0e-5x₃²x₄⁶ + 6.0e-5x₃⁴x₄⁴ + 4.0e-5x₃⁶x₄² + 1.0e-5x₃⁸ + 4.0e-5x₂²x₄⁶ + 0.00012x₂²x₃²x₄⁴ + 0.00012x₂²x₃⁴x₄² + 4.0e-5x₂²x₃⁶ + 6.0e-5x₂⁴x₄⁴ + 0.00012x₂⁴x₃²x₄² + 6.0e-5x₂⁴x₃⁴ + 4.0e-5x₂⁶x₄² + 4.0e-5x₂⁶x₃² + 1.0e-5x₂⁸ + 4.0e-5x₁²x₄⁶ + 0.00012x₁²x₃²x₄⁴ + 0.00012x₁²x₃⁴x₄² + 4.0e-5x₁²x₃⁶ + 0.00012x₁²x₂²x₄⁴ + 0.00024x₁²x₂²x₃²x₄² + 0.00012x₁²x₂²x₃⁴ + 0.00012x₁²x₂⁴x₄² + 0.00012x₁²x₂⁴x₃² + 4.0e-5x₁²x₂⁶ + 6.0e-5x₁⁴x₄⁴ + 0.00012x₁⁴x₃²x₄² + 6.0e-5x₁⁴x₃⁴ + 0.00012x₁⁴x₂²x₄² + 0.00012x₁⁴x₂²x₃² + 6.0e-5x₁⁴x₂⁴ + 4.0e-5x₁⁶x₄² + 4.0e-5x₁⁶x₃² + 4.0e-5x₁⁶x₂² + 1.0e-5x₁⁸
Objective was scaled by the prefactor 1.0 + x₄² + x₃² + x₂² + x₁²
3 equality constraints
1: x₁ = 0
2: -x₃ + x₂ = 0
3: -x₄ + x₃ = 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = 5e-5
            end
        end
    end
end

@testset "POP 16 (Example A.4)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem((x[1] +1)^2 + x[2]^2, zero=[x[1]^3-x[2]^2], noncompact=(1e-5, 5))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.00001 + 2.0x₁ + 6.00007x₂² + 6.00007x₁² + 10.0x₁x₂² + 10.0x₁³ + 15.00021x₂⁴ + 30.00042x₁²x₂² + 15.00021x₁⁴ + 20.0x₁x₂⁴ + 40.0x₁³x₂² + 20.0x₁⁵ + 20.00035x₂⁶ + 60.00105x₁²x₂⁴ + 60.00105x₁⁴x₂² + 20.00035x₁⁶ + 20.0x₁x₂⁶ + 60.0x₁³x₂⁴ + 60.0x₁⁵x₂² + 20.0x₁⁷ + 15.00035x₂⁸ + 60.0014x₁²x₂⁶ + 90.0021x₁⁴x₂⁴ + 60.0014x₁⁶x₂² + 15.00035x₁⁸ + 10.0x₁x₂⁸ + 40.0x₁³x₂⁶ + 60.0x₁⁵x₂⁴ + 40.0x₁⁷x₂² + 10.0x₁⁹ + 6.00021x₂¹⁰ + 30.00105x₁²x₂⁸ + 60.0021x₁⁴x₂⁶ + 60.0021x₁⁶x₂⁴ + 30.00105x₁⁸x₂² + 6.00021x₁¹⁰ + 2.0x₁x₂¹⁰ + 10.0x₁³x₂⁸ + 20.0x₁⁵x₂⁶ + 20.0x₁⁷x₂⁴ + 10.0x₁⁹x₂² + 2.0x₁¹¹ + 1.00007x₂¹² + 6.00042x₁²x₂¹⁰ + 15.00105x₁⁴x₂⁸ + 20.0014x₁⁶x₂⁶ + 15.00105x₁⁸x₂⁴ + 6.00042x₁¹⁰x₂² + 1.00007x₁¹² + 1.0e-5x₂¹⁴ + 7.0e-5x₁²x₂¹² + 0.00021x₁⁴x₂¹⁰ + 0.00035x₁⁶x₂⁸ + 0.00035x₁⁸x₂⁶ + 0.00021x₁¹⁰x₂⁴ + 7.0e-5x₁¹²x₂² + 1.0e-5x₁¹⁴
Objective was scaled by the prefactor 1.0 + 5.0x₂² + 5.0x₁² + 10.0x₂⁴ + 20.0x₁²x₂² + 10.0x₁⁴ + 10.0x₂⁶ + 30.0x₁²x₂⁴ + 30.0x₁⁴x₂² + 10.0x₁⁶ + 5.0x₂⁸ + 20.0x₁²x₂⁶ + 30.0x₁⁴x₂⁴ + 20.0x₁⁶x₂² + 5.0x₁⁸ + x₂¹⁰ + 5.0x₁²x₂⁸ + 10.0x₁⁴x₂⁶ + 10.0x₁⁶x₂⁴ + 5.0x₁⁸x₂² + x₁¹⁰
1 equality constraint
1: -x₂² + x₁³ = 0"
    # this problem is very ill-defined, none of the solver will report success

    Random.seed!(1234)
    prob = poly_problem((x[1] +1)^2 + x[2]^2, zero=[x[1]^3-x[2]^2], noncompact=(1e-5, 10))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.00001 + 2.0x₁ + 11.00012x₂² + 11.00012x₁² + 20.0x₁x₂² + 20.0x₁³ + 55.00066x₂⁴ + 110.00132x₁²x₂² + 55.00066x₁⁴ + 90.0x₁x₂⁴ + 180.0x₁³x₂² + 90.0x₁⁵ + 165.0022x₂⁶ + 495.0066x₁²x₂⁴ + 495.0066x₁⁴x₂² + 165.0022x₁⁶ + 240.0x₁x₂⁶ + 720.0x₁³x₂⁴ + 720.0x₁⁵x₂² + 240.0x₁⁷ + 330.00495x₂⁸ + 1320.0198x₁²x₂⁶ + 1980.0297x₁⁴x₂⁴ + 1320.0198x₁⁶x₂² + 330.00495x₁⁸ + 420.0x₁x₂⁸ + 1680.0x₁³x₂⁶ + 2520.0x₁⁵x₂⁴ + 1680.0x₁⁷x₂² + 420.0x₁⁹ + 462.00792x₂¹⁰ + 2310.0396x₁²x₂⁸ + 4620.0792x₁⁴x₂⁶ + 4620.0792x₁⁶x₂⁴ + 2310.0396x₁⁸x₂² + 462.00792x₁¹⁰ + 504.0x₁x₂¹⁰ + 2520.0x₁³x₂⁸ + 5040.0x₁⁵x₂⁶ + 5040.0x₁⁷x₂⁴ + 2520.0x₁⁹x₂² + 504.0x₁¹¹ + 462.00924x₂¹² + 2772.05544x₁²x₂¹⁰ + 6930.1386x₁⁴x₂⁸ + 9240.1848x₁⁶x₂⁶ + 6930.1386x₁⁸x₂⁴ + 2772.05544x₁¹⁰x₂² + 462.00924x₁¹² + 420.0x₁x₂¹² + 2520.0x₁³x₂¹⁰ + 6300.0x₁⁵x₂⁸ + 8400.0x₁⁷x₂⁶ + 6300.0x₁⁹x₂⁴ + 2520.0x₁¹¹x₂² + 420.0x₁¹³ + 330.00792x₂¹⁴ + 2310.05544x₁²x₂¹² + 6930.16632x₁⁴x₂¹⁰ + 11550.2772x₁⁶x₂⁸ + 11550.2772x₁⁸x₂⁶ + 6930.16632x₁¹⁰x₂⁴ + 2310.05544x₁¹²x₂² + 330.00792x₁¹⁴ + 240.0x₁x₂¹⁴ + 1680.0x₁³x₂¹² + 5040.0x₁⁵x₂¹⁰ + 8400.0x₁⁷x₂⁸ + 8400.0x₁⁹x₂⁶ + 5040.0x₁¹¹x₂⁴ + 1680.0x₁¹³x₂² + 240.0x₁¹⁵ + 165.00495x₂¹⁶ + 1320.0396x₁²x₂¹⁴ + 4620.1386x₁⁴x₂¹² + 9240.2772x₁⁶x₂¹⁰ + 11550.3465x₁⁸x₂⁸ + 9240.2772x₁¹⁰x₂⁶ + 4620.1386x₁¹²x₂⁴ + 1320.0396x₁¹⁴x₂² + 165.00495x₁¹⁶ + 90.0x₁x₂¹⁶ + 720.0x₁³x₂¹⁴ + 2520.0x₁⁵x₂¹² + 5040.0x₁⁷x₂¹⁰ + 6300.0x₁⁹x₂⁸ + 5040.0x₁¹¹x₂⁶ + 2520.0x₁¹³x₂⁴ + 720.0x₁¹⁵x₂² + 90.0x₁¹⁷ + 55.0022x₂¹⁸ + 495.0198x₁²x₂¹⁶ + 1980.0792x₁⁴x₂¹⁴ + 4620.1848x₁⁶x₂¹² + 6930.2772x₁⁸x₂¹⁰ + 6930.2772x₁¹⁰x₂⁸ + 4620.1848x₁¹²x₂⁶ + 1980.0792x₁¹⁴x₂⁴ + 495.0198x₁¹⁶x₂² + 55.0022x₁¹⁸ + 20.0x₁x₂¹⁸ + 180.0x₁³x₂¹⁶ + 720.0x₁⁵x₂¹⁴ + 1680.0x₁⁷x₂¹² + 2520.0x₁⁹x₂¹⁰ + 2520.0x₁¹¹x₂⁸ + 1680.0x₁¹³x₂⁶ + 720.0x₁¹⁵x₂⁴ + 180.0x₁¹⁷x₂² + 20.0x₁¹⁹ + 11.00066x₂²⁰ + 110.0066x₁²x₂¹⁸ + 495.0297x₁⁴x₂¹⁶ + 1320.0792x₁⁶x₂¹⁴ + 2310.1386x₁⁸x₂¹² + 2772.16632x₁¹⁰x₂¹⁰ + 2310.1386x₁¹²x₂⁸ + 1320.0792x₁¹⁴x₂⁶ + 495.0297x₁¹⁶x₂⁴ + 110.0066x₁¹⁸x₂² + 11.00066x₁²⁰ + 2.0x₁x₂²⁰ + 20.0x₁³x₂¹⁸ + 90.0x₁⁵x₂¹⁶ + 240.0x₁⁷x₂¹⁴ + 420.0x₁⁹x₂¹² + 504.0x₁¹¹x₂¹⁰ + 420.0x₁¹³x₂⁸ + 240.0x₁¹⁵x₂⁶ + 90.0x₁¹⁷x₂⁴ + 20.0x₁¹⁹x₂² + 2.0x₁²¹ + 1.00012x₂²² + 11.00132x₁²x₂²⁰ + 55.0066x₁⁴x₂¹⁸ + 165.0198x₁⁶x₂¹⁶ + 330.0396x₁⁸x₂¹⁴ + 462.05544x₁¹⁰x₂¹² + 462.05544x₁¹²x₂¹⁰ + 330.0396x₁¹⁴x₂⁸ + 165.0198x₁¹⁶x₂⁶ + 55.0066x₁¹⁸x₂⁴ + 11.00132x₁²⁰x₂² + 1.00012x₁²² + 1.0e-5x₂²⁴ + 0.00012x₁²x₂²² + 0.00066x₁⁴x₂²⁰ + 0.0022x₁⁶x₂¹⁸ + 0.00495x₁⁸x₂¹⁶ + 0.00792x₁¹⁰x₂¹⁴ + 0.00924x₁¹²x₂¹² + 0.00792x₁¹⁴x₂¹⁰ + 0.00495x₁¹⁶x₂⁸ + 0.0022x₁¹⁸x₂⁶ + 0.00066x₁²⁰x₂⁴ + 0.00012x₁²²x₂² + 1.0e-5x₁²⁴
Objective was scaled by the prefactor 1.0 + 10.0x₂² + 10.0x₁² + 45.0x₂⁴ + 90.0x₁²x₂² + 45.0x₁⁴ + 120.0x₂⁶ + 360.0x₁²x₂⁴ + 360.0x₁⁴x₂² + 120.0x₁⁶ + 210.0x₂⁸ + 840.0x₁²x₂⁶ + 1260.0x₁⁴x₂⁴ + 840.0x₁⁶x₂² + 210.0x₁⁸ + 252.0x₂¹⁰ + 1260.0x₁²x₂⁸ + 2520.0x₁⁴x₂⁶ + 2520.0x₁⁶x₂⁴ + 1260.0x₁⁸x₂² + 252.0x₁¹⁰ + 210.0x₂¹² + 1260.0x₁²x₂¹⁰ + 3150.0x₁⁴x₂⁸ + 4200.0x₁⁶x₂⁶ + 3150.0x₁⁸x₂⁴ + 1260.0x₁¹⁰x₂² + 210.0x₁¹² + 120.0x₂¹⁴ + 840.0x₁²x₂¹² + 2520.0x₁⁴x₂¹⁰ + 4200.0x₁⁶x₂⁸ + 4200.0x₁⁸x₂⁶ + 2520.0x₁¹⁰x₂⁴ + 840.0x₁¹²x₂² + 120.0x₁¹⁴ + 45.0x₂¹⁶ + 360.0x₁²x₂¹⁴ + 1260.0x₁⁴x₂¹² + 2520.0x₁⁶x₂¹⁰ + 3150.0x₁⁸x₂⁸ + 2520.0x₁¹⁰x₂⁶ + 1260.0x₁¹²x₂⁴ + 360.0x₁¹⁴x₂² + 45.0x₁¹⁶ + 10.0x₂¹⁸ + 90.0x₁²x₂¹⁶ + 360.0x₁⁴x₂¹⁴ + 840.0x₁⁶x₂¹² + 1260.0x₁⁸x₂¹⁰ + 1260.0x₁¹⁰x₂⁸ + 840.0x₁¹²x₂⁶ + 360.0x₁¹⁴x₂⁴ + 90.0x₁¹⁶x₂² + 10.0x₁¹⁸ + x₂²⁰ + 10.0x₁²x₂¹⁸ + 45.0x₁⁴x₂¹⁶ + 120.0x₁⁶x₂¹⁴ + 210.0x₁⁸x₂¹² + 252.0x₁¹⁰x₂¹⁰ + 210.0x₁¹²x₂⁸ + 120.0x₁¹⁴x₂⁶ + 45.0x₁¹⁶x₂⁴ + 10.0x₁¹⁸x₂² + x₁²⁰
1 equality constraint
1: -x₂² + x₁³ = 0"
    # same problem

    Random.seed!(1234)
    prob = poly_problem((x[1] +1)^2 + x[2]^2, zero=[x[1]^3-x[2]^2], noncompact=(1e-5, 15))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.00001 + 2.0x₁ + 16.00017x₂² + 16.00017x₁² + 30.0x₁x₂² + 30.0x₁³ + 120.00136x₂⁴ + 240.00272x₁²x₂² + 120.00136x₁⁴ + 210.0x₁x₂⁴ + 420.0x₁³x₂² + 210.0x₁⁵ + 560.0068x₂⁶ + 1680.0204x₁²x₂⁴ + 1680.0204x₁⁴x₂² + 560.0068x₁⁶ + 910.0x₁x₂⁶ + 2730.0x₁³x₂⁴ + 2730.0x₁⁵x₂² + 910.0x₁⁷ + 1820.0238x₂⁸ + 7280.0952x₁²x₂⁶ + 10920.1428x₁⁴x₂⁴ + 7280.0952x₁⁶x₂² + 1820.0238x₁⁸ + 2730.0x₁x₂⁸ + 10920.0x₁³x₂⁶ + 16380.0x₁⁵x₂⁴ + 10920.0x₁⁷x₂² + 2730.0x₁⁹ + 4368.06188x₂¹⁰ + 21840.3094x₁²x₂⁸ + 43680.6188x₁⁴x₂⁶ + 43680.6188x₁⁶x₂⁴ + 21840.3094x₁⁸x₂² + 4368.06188x₁¹⁰ + 6006.0x₁x₂¹⁰ + 30030.0x₁³x₂⁸ + 60060.0x₁⁵x₂⁶ + 60060.0x₁⁷x₂⁴ + 30030.0x₁⁹x₂² + 6006.0x₁¹¹ + 8008.12376x₂¹² + 48048.74256x₁²x₂¹⁰ + 120121.8564x₁⁴x₂⁸ + 160162.4752x₁⁶x₂⁶ + 120121.8564x₁⁸x₂⁴ + 48048.74256x₁¹⁰x₂² + 8008.12376x₁¹² + 10010.0x₁x₂¹² + 60060.0x₁³x₂¹⁰ + 150150.0x₁⁵x₂⁸ + 200200.0x₁⁷x₂⁶ + 150150.0x₁⁹x₂⁴ + 60060.0x₁¹¹x₂² + 10010.0x₁¹³ + 11440.19448x₂¹⁴ + 80081.36136x₁²x₂¹² + 240244.08408x₁⁴x₂¹⁰ + 400406.8068x₁⁶x₂⁸ + 400406.8068x₁⁸x₂⁶ + 240244.08408x₁¹⁰x₂⁴ + 80081.36136x₁¹²x₂² + 11440.19448x₁¹⁴ + 12870.0x₁x₂¹⁴ + 90090.0x₁³x₂¹² + 270270.0x₁⁵x₂¹⁰ + 450450.0x₁⁷x₂⁸ + 450450.0x₁⁹x₂⁶ + 270270.0x₁¹¹x₂⁴ + 90090.0x₁¹³x₂² + 12870.0x₁¹⁵ + 12870.2431x₂¹⁶ + 102961.9448x₁²x₂¹⁴ + 360366.8068x₁⁴x₂¹² + 720733.6136x₁⁶x₂¹⁰ + 900917.017x₁⁸x₂⁸ + 720733.6136x₁¹⁰x₂⁶ + 360366.8068x₁¹²x₂⁴ + 102961.9448x₁¹⁴x₂² + 12870.2431x₁¹⁶ + 12870.0x₁x₂¹⁶ + 102960.0x₁³x₂¹⁴ + 360360.0x₁⁵x₂¹² + 720720.0x₁⁷x₂¹⁰ + 900900.0x₁⁹x₂⁸ + 720720.0x₁¹¹x₂⁶ + 360360.0x₁¹³x₂⁴ + 102960.0x₁¹⁵x₂² + 12870.0x₁¹⁷ + 11440.2431x₂¹⁸ + 102962.1879x₁²x₂¹⁶ + 411848.7516x₁⁴x₂¹⁴ + 960980.4204x₁⁶x₂¹² + 1.4414706306e6x₁⁸x₂¹⁰ + 1.4414706306e6x₁¹⁰x₂⁸ + 960980.4204x₁¹²x₂⁶ + 411848.7516x₁¹⁴x₂⁴ + 102962.1879x₁¹⁶x₂² + 11440.2431x₁¹⁸ + 10010.0x₁x₂¹⁸ + 90090.0x₁³x₂¹⁶ + 360360.0x₁⁵x₂¹⁴ + 840840.0x₁⁷x₂¹² + 1.26126e6x₁⁹x₂¹⁰ + 1.26126e6x₁¹¹x₂⁸ + 840840.0x₁¹³x₂⁶ + 360360.0x₁¹⁵x₂⁴ + 90090.0x₁¹⁷x₂² + 10010.0x₁¹⁹ + 8008.19448x₂²⁰ + 80081.9448x₁²x₂¹⁸ + 360368.7516x₁⁴x₂¹⁶ + 960983.3376x₁⁶x₂¹⁴ + 1.6817208408e6x₁⁸x₂¹² + 2.01806500896e6x₁¹⁰x₂¹⁰ + 1.6817208408e6x₁¹²x₂⁸ + 960983.3376x₁¹⁴x₂⁶ + 360368.7516x₁¹⁶x₂⁴ + 80081.9448x₁¹⁸x₂² + 8008.19448x₁²⁰ + 6006.0x₁x₂²⁰ + 60060.0x₁³x₂¹⁸ + 270270.0x₁⁵x₂¹⁶ + 720720.0x₁⁷x₂¹⁴ + 1.26126e6x₁⁹x₂¹² + 1.513512e6x₁¹¹x₂¹⁰ + 1.26126e6x₁¹³x₂⁸ + 720720.0x₁¹⁵x₂⁶ + 270270.0x₁¹⁷x₂⁴ + 60060.0x₁¹⁹x₂² + 6006.0x₁²¹ + 4368.12376x₂²² + 48049.36136x₁²x₂²⁰ + 240246.8068x₁⁴x₂¹⁸ + 720740.4204x₁⁶x₂¹⁶ + 1.4414808408e6x₁⁸x₂¹⁴ + 2.01807317712e6x₁¹⁰x₂¹² + 2.01807317712e6x₁¹²x₂¹⁰ + 1.4414808408e6x₁¹⁴x₂⁸ + 720740.4204x₁¹⁶x₂⁶ + 240246.8068x₁¹⁸x₂⁴ + 48049.36136x₁²⁰x₂² + 4368.12376x₁²² + 2730.0x₁x₂²² + 30030.0x₁³x₂²⁰ + 150150.0x₁⁵x₂¹⁸ + 450450.0x₁⁷x₂¹⁶ + 900900.0x₁⁹x₂¹⁴ + 1.26126e6x₁¹¹x₂¹² + 1.26126e6x₁¹³x₂¹⁰ + 900900.0x₁¹⁵x₂⁸ + 450450.0x₁¹⁷x₂⁶ + 150150.0x₁¹⁹x₂⁴ + 30030.0x₁²¹x₂² + 2730.0x₁²³ + 1820.06188x₂²⁴ + 21840.74256x₁²x₂²² + 120124.08408x₁⁴x₂²⁰ + 400413.6136x₁⁶x₂¹⁸ + 900930.6306x₁⁸x₂¹⁶ + 1.44148900896e6x₁¹⁰x₂¹⁴ + 1.68173717712e6x₁¹²x₂¹² + 1.44148900896e6x₁¹⁴x₂¹⁰ + 900930.6306x₁¹⁶x₂⁸ + 400413.6136x₁¹⁸x₂⁶ + 120124.08408x₁²⁰x₂⁴ + 21840.74256x₁²²x₂² + 1820.06188x₁²⁴ + 910.0x₁x₂²⁴ + 10920.0x₁³x₂²² + 60060.0x₁⁵x₂²⁰ + 200200.0x₁⁷x₂¹⁸ + 450450.0x₁⁹x₂¹⁶ + 720720.0x₁¹¹x₂¹⁴ + 840840.0x₁¹³x₂¹² + 720720.0x₁¹⁵x₂¹⁰ + 450450.0x₁¹⁷x₂⁸ + 200200.0x₁¹⁹x₂⁶ + 60060.0x₁²¹x₂⁴ + 10920.0x₁²³x₂² + 910.0x₁²⁵ + 560.0238x₂²⁶ + 7280.3094x₁²x₂²⁴ + 43681.8564x₁⁴x₂²² + 160166.8068x₁⁶x₂²⁰ + 400417.017x₁⁸x₂¹⁸ + 720750.6306x₁¹⁰x₂¹⁶ + 961000.8408x₁¹²x₂¹⁴ + 961000.8408x₁¹⁴x₂¹² + 720750.6306x₁¹⁶x₂¹⁰ + 400417.017x₁¹⁸x₂⁸ + 160166.8068x₁²⁰x₂⁶ + 43681.8564x₁²²x₂⁴ + 7280.3094x₁²⁴x₂² + 560.0238x₁²⁶ + 210.0x₁x₂²⁶ + 2730.0x₁³x₂²⁴ + 16380.0x₁⁵x₂²² + 60060.0x₁⁷x₂²⁰ + 150150.0x₁⁹x₂¹⁸ + 270270.0x₁¹¹x₂¹⁶ + 360360.0x₁¹³x₂¹⁴ + 360360.0x₁¹⁵x₂¹² + 270270.0x₁¹⁷x₂¹⁰ + 150150.0x₁¹⁹x₂⁸ + 60060.0x₁²¹x₂⁶ + 16380.0x₁²³x₂⁴ + 2730.0x₁²⁵x₂² + 210.0x₁²⁷ + 120.0068x₂²⁸ + 1680.0952x₁²x₂²⁶ + 10920.6188x₁⁴x₂²⁴ + 43682.4752x₁⁶x₂²² + 120126.8068x₁⁸x₂²⁰ + 240253.6136x₁¹⁰x₂¹⁸ + 360380.4204x₁¹²x₂¹⁶ + 411863.3376x₁¹⁴x₂¹⁴ + 360380.4204x₁¹⁶x₂¹² + 240253.6136x₁¹⁸x₂¹⁰ + 120126.8068x₁²⁰x₂⁸ + 43682.4752x₁²²x₂⁶ + 10920.6188x₁²⁴x₂⁴ + 1680.0952x₁²⁶x₂² + 120.0068x₁²⁸ + 30.0x₁x₂²⁸ + 420.0x₁³x₂²⁶ + 2730.0x₁⁵x₂²⁴ + 10920.0x₁⁷x₂²² + 30030.0x₁⁹x₂²⁰ + 60060.0x₁¹¹x₂¹⁸ + 90090.0x₁¹³x₂¹⁶ + 102960.0x₁¹⁵x₂¹⁴ + 90090.0x₁¹⁷x₂¹² + 60060.0x₁¹⁹x₂¹⁰ + 30030.0x₁²¹x₂⁸ + 10920.0x₁²³x₂⁶ + 2730.0x₁²⁵x₂⁴ + 420.0x₁²⁷x₂² + 30.0x₁²⁹ + 16.00136x₂³⁰ + 240.0204x₁²x₂²⁸ + 1680.1428x₁⁴x₂²⁶ + 7280.6188x₁⁶x₂²⁴ + 21841.8564x₁⁸x₂²² + 48052.08408x₁¹⁰x₂²⁰ + 80086.8068x₁¹²x₂¹⁸ + 102968.7516x₁¹⁴x₂¹⁶ + 102968.7516x₁¹⁶x₂¹⁴ + 80086.8068x₁¹⁸x₂¹² + 48052.08408x₁²⁰x₂¹⁰ + 21841.8564x₁²²x₂⁸ + 7280.6188x₁²⁴x₂⁶ + 1680.1428x₁²⁶x₂⁴ + 240.0204x₁²⁸x₂² + 16.00136x₁³⁰ + 2.0x₁x₂³⁰ + 30.0x₁³x₂²⁸ + 210.0x₁⁵x₂²⁶ + 910.0x₁⁷x₂²⁴ + 2730.0x₁⁹x₂²² + 6006.0x₁¹¹x₂²⁰ + 10010.0x₁¹³x₂¹⁸ + 12870.0x₁¹⁵x₂¹⁶ + 12870.0x₁¹⁷x₂¹⁴ + 10010.0x₁¹⁹x₂¹² + 6006.0x₁²¹x₂¹⁰ + 2730.0x₁²³x₂⁸ + 910.0x₁²⁵x₂⁶ + 210.0x₁²⁷x₂⁴ + 30.0x₁²⁹x₂² + 2.0x₁³¹ + 1.00017x₂³² + 16.00272x₁²x₂³⁰ + 120.0204x₁⁴x₂²⁸ + 560.0952x₁⁶x₂²⁶ + 1820.3094x₁⁸x₂²⁴ + 4368.74256x₁¹⁰x₂²² + 8009.36136x₁¹²x₂²⁰ + 11441.9448x₁¹⁴x₂¹⁸ + 12872.1879x₁¹⁶x₂¹⁶ + 11441.9448x₁¹⁸x₂¹⁴ + 8009.36136x₁²⁰x₂¹² + 4368.74256x₁²²x₂¹⁰ + 1820.3094x₁²⁴x₂⁸ + 560.0952x₁²⁶x₂⁶ + 120.0204x₁²⁸x₂⁴ + 16.00272x₁³⁰x₂² + 1.00017x₁³² + 1.0e-5x₂³⁴ + 0.00017x₁²x₂³² + 0.00136x₁⁴x₂³⁰ + 0.0068x₁⁶x₂²⁸ + 0.0238x₁⁸x₂²⁶ + 0.06188x₁¹⁰x₂²⁴ + 0.12376x₁¹²x₂²² + 0.19448x₁¹⁴x₂²⁰ + 0.2431x₁¹⁶x₂¹⁸ + 0.2431x₁¹⁸x₂¹⁶ + 0.19448x₁²⁰x₂¹⁴ + 0.12376x₁²²x₂¹² + 0.06188x₁²⁴x₂¹⁰ + 0.0238x₁²⁶x₂⁸ + 0.0068x₁²⁸x₂⁶ + 0.00136x₁³⁰x₂⁴ + 0.00017x₁³²x₂² + 1.0e-5x₁³⁴
Objective was scaled by the prefactor 1.0 + 15.0x₂² + 15.0x₁² + 105.0x₂⁴ + 210.0x₁²x₂² + 105.0x₁⁴ + 455.0x₂⁶ + 1365.0x₁²x₂⁴ + 1365.0x₁⁴x₂² + 455.0x₁⁶ + 1365.0x₂⁸ + 5460.0x₁²x₂⁶ + 8190.0x₁⁴x₂⁴ + 5460.0x₁⁶x₂² + 1365.0x₁⁸ + 3003.0x₂¹⁰ + 15015.0x₁²x₂⁸ + 30030.0x₁⁴x₂⁶ + 30030.0x₁⁶x₂⁴ + 15015.0x₁⁸x₂² + 3003.0x₁¹⁰ + 5005.0x₂¹² + 30030.0x₁²x₂¹⁰ + 75075.0x₁⁴x₂⁸ + 100100.0x₁⁶x₂⁶ + 75075.0x₁⁸x₂⁴ + 30030.0x₁¹⁰x₂² + 5005.0x₁¹² + 6435.0x₂¹⁴ + 45045.0x₁²x₂¹² + 135135.0x₁⁴x₂¹⁰ + 225225.0x₁⁶x₂⁸ + 225225.0x₁⁸x₂⁶ + 135135.0x₁¹⁰x₂⁴ + 45045.0x₁¹²x₂² + 6435.0x₁¹⁴ + 6435.0x₂¹⁶ + 51480.0x₁²x₂¹⁴ + 180180.0x₁⁴x₂¹² + 360360.0x₁⁶x₂¹⁰ + 450450.0x₁⁸x₂⁸ + 360360.0x₁¹⁰x₂⁶ + 180180.0x₁¹²x₂⁴ + 51480.0x₁¹⁴x₂² + 6435.0x₁¹⁶ + 5005.0x₂¹⁸ + 45045.0x₁²x₂¹⁶ + 180180.0x₁⁴x₂¹⁴ + 420420.0x₁⁶x₂¹² + 630630.0x₁⁸x₂¹⁰ + 630630.0x₁¹⁰x₂⁸ + 420420.0x₁¹²x₂⁶ + 180180.0x₁¹⁴x₂⁴ + 45045.0x₁¹⁶x₂² + 5005.0x₁¹⁸ + 3003.0x₂²⁰ + 30030.0x₁²x₂¹⁸ + 135135.0x₁⁴x₂¹⁶ + 360360.0x₁⁶x₂¹⁴ + 630630.0x₁⁸x₂¹² + 756756.0x₁¹⁰x₂¹⁰ + 630630.0x₁¹²x₂⁸ + 360360.0x₁¹⁴x₂⁶ + 135135.0x₁¹⁶x₂⁴ + 30030.0x₁¹⁸x₂² + 3003.0x₁²⁰ + 1365.0x₂²² + 15015.0x₁²x₂²⁰ + 75075.0x₁⁴x₂¹⁸ + 225225.0x₁⁶x₂¹⁶ + 450450.0x₁⁸x₂¹⁴ + 630630.0x₁¹⁰x₂¹² + 630630.0x₁¹²x₂¹⁰ + 450450.0x₁¹⁴x₂⁸ + 225225.0x₁¹⁶x₂⁶ + 75075.0x₁¹⁸x₂⁴ + 15015.0x₁²⁰x₂² + 1365.0x₁²² + 455.0x₂²⁴ + 5460.0x₁²x₂²² + 30030.0x₁⁴x₂²⁰ + 100100.0x₁⁶x₂¹⁸ + 225225.0x₁⁸x₂¹⁶ + 360360.0x₁¹⁰x₂¹⁴ + 420420.0x₁¹²x₂¹² + 360360.0x₁¹⁴x₂¹⁰ + 225225.0x₁¹⁶x₂⁸ + 100100.0x₁¹⁸x₂⁶ + 30030.0x₁²⁰x₂⁴ + 5460.0x₁²²x₂² + 455.0x₁²⁴ + 105.0x₂²⁶ + 1365.0x₁²x₂²⁴ + 8190.0x₁⁴x₂²² + 30030.0x₁⁶x₂²⁰ + 75075.0x₁⁸x₂¹⁸ + 135135.0x₁¹⁰x₂¹⁶ + 180180.0x₁¹²x₂¹⁴ + 180180.0x₁¹⁴x₂¹² + 135135.0x₁¹⁶x₂¹⁰ + 75075.0x₁¹⁸x₂⁸ + 30030.0x₁²⁰x₂⁶ + 8190.0x₁²²x₂⁴ + 1365.0x₁²⁴x₂² + 105.0x₁²⁶ + 15.0x₂²⁸ + 210.0x₁²x₂²⁶ + 1365.0x₁⁴x₂²⁴ + 5460.0x₁⁶x₂²² + 15015.0x₁⁸x₂²⁰ + 30030.0x₁¹⁰x₂¹⁸ + 45045.0x₁¹²x₂¹⁶ + 51480.0x₁¹⁴x₂¹⁴ + 45045.0x₁¹⁶x₂¹² + 30030.0x₁¹⁸x₂¹⁰ + 15015.0x₁²⁰x₂⁸ + 5460.0x₁²²x₂⁶ + 1365.0x₁²⁴x₂⁴ + 210.0x₁²⁶x₂² + 15.0x₁²⁸ + x₂³⁰ + 15.0x₁²x₂²⁸ + 105.0x₁⁴x₂²⁶ + 455.0x₁⁶x₂²⁴ + 1365.0x₁⁸x₂²² + 3003.0x₁¹⁰x₂²⁰ + 5005.0x₁¹²x₂¹⁸ + 6435.0x₁¹⁴x₂¹⁶ + 6435.0x₁¹⁶x₂¹⁴ + 5005.0x₁¹⁸x₂¹² + 3003.0x₁²⁰x₂¹⁰ + 1365.0x₁²²x₂⁸ + 455.0x₁²⁴x₂⁶ + 105.0x₁²⁶x₂⁴ + 15.0x₁²⁸x₂² + x₁³⁰
1 equality constraint
1: -x₂² + x₁³ = 0"
    if optimize
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, prob).objective ≈ 0.9858 atol = 1e-4
        :HypatiaMoment ∈ solvers && @test poly_optimize(:HypatiaMoment, prob).objective ≈ 0.9858 atol = 1e-4
        # Hypatia is the only one of the others that gives something close to a good result (30s)
    end

    Random.seed!(1234)
    prob = poly_problem((x[1] +1)^2 + x[2]^2, zero=[x[1]^3-x[2]^2], noncompact=(1e-5, 20))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.00001 + 2.0x₁ + 21.00022x₂² + 21.00022x₁² + 40.0x₁x₂² + 40.0x₁³ + 210.00231x₂⁴ + 420.00462x₁²x₂² + 210.00231x₁⁴ + 380.0x₁x₂⁴ + 760.0x₁³x₂² + 380.0x₁⁵ + 1330.0154x₂⁶ + 3990.0462x₁²x₂⁴ + 3990.0462x₁⁴x₂² + 1330.0154x₁⁶ + 2280.0x₁x₂⁶ + 6840.0x₁³x₂⁴ + 6840.0x₁⁵x₂² + 2280.0x₁⁷ + 5985.07315x₂⁸ + 23940.2926x₁²x₂⁶ + 35910.4389x₁⁴x₂⁴ + 23940.2926x₁⁶x₂² + 5985.07315x₁⁸ + 9690.0x₁x₂⁸ + 38760.0x₁³x₂⁶ + 58140.0x₁⁵x₂⁴ + 38760.0x₁⁷x₂² + 9690.0x₁⁹ + 20349.26334x₂¹⁰ + 101746.3167x₁²x₂⁸ + 203492.6334x₁⁴x₂⁶ + 203492.6334x₁⁶x₂⁴ + 101746.3167x₁⁸x₂² + 20349.26334x₁¹⁰ + 31008.0x₁x₂¹⁰ + 155040.0x₁³x₂⁸ + 310080.0x₁⁵x₂⁶ + 310080.0x₁⁷x₂⁴ + 155040.0x₁⁹x₂² + 31008.0x₁¹¹ + 54264.74613x₂¹² + 325588.47678x₁²x₂¹⁰ + 813971.19195x₁⁴x₂⁸ + 1.0852949226e6x₁⁶x₂⁶ + 813971.19195x₁⁸x₂⁴ + 325588.47678x₁¹⁰x₂² + 54264.74613x₁¹² + 77520.0x₁x₂¹² + 465120.0x₁³x₂¹⁰ + 1.1628e6x₁⁵x₂⁸ + 1.5504e6x₁⁷x₂⁶ + 1.1628e6x₁⁹x₂⁴ + 465120.0x₁¹¹x₂² + 77520.0x₁¹³ + 116281.70544x₂¹⁴ + 813971.93808x₁²x₂¹² + 2.44191581424e6x₁⁴x₂¹⁰ + 4.0698596904e6x₁⁶x₂⁸ + 4.0698596904e6x₁⁸x₂⁶ + 2.44191581424e6x₁¹⁰x₂⁴ + 813971.93808x₁¹²x₂² + 116281.70544x₁¹⁴ + 155040.0x₁x₂¹⁴ + 1.08528e6x₁³x₂¹² + 3.25584e6x₁⁵x₂¹⁰ + 5.4264e6x₁⁷x₂⁸ + 5.4264e6x₁⁹x₂⁶ + 3.25584e6x₁¹¹x₂⁴ + 1.08528e6x₁¹³x₂² + 155040.0x₁¹⁵ + 203493.1977x₂¹⁶ + 1.6279455816e6x₁²x₂¹⁴ + 5.6978095356e6x₁⁴x₂¹² + 1.13956190712e7x₁⁶x₂¹⁰ + 1.4244523839e7x₁⁸x₂⁸ + 1.13956190712e7x₁¹⁰x₂⁶ + 5.6978095356e6x₁¹²x₂⁴ + 1.6279455816e6x₁¹⁴x₂² + 203493.1977x₁¹⁶ + 251940.0x₁x₂¹⁶ + 2.01552e6x₁³x₂¹⁴ + 7.05432e6x₁⁵x₂¹² + 1.410864e7x₁⁷x₂¹⁰ + 1.76358e7x₁⁹x₂⁸ + 1.410864e7x₁¹¹x₂⁶ + 7.05432e6x₁¹³x₂⁴ + 2.01552e6x₁¹⁵x₂² + 251940.0x₁¹⁷ + 293934.9742x₂¹⁸ + 2.6454147678e6x₁²x₂¹⁶ + 1.05816590712e7x₁⁴x₂¹⁴ + 2.46905378328e7x₁⁶x₂¹² + 3.70358067492e7x₁⁸x₂¹⁰ + 3.70358067492e7x₁¹⁰x₂⁸ + 2.46905378328e7x₁¹²x₂⁶ + 1.05816590712e7x₁¹⁴x₂⁴ + 2.6454147678e6x₁¹⁶x₂² + 293934.9742x₁¹⁸ + 335920.0x₁x₂¹⁸ + 3.02328e6x₁³x₂¹⁶ + 1.209312e7x₁⁵x₂¹⁴ + 2.821728e7x₁⁷x₂¹² + 4.232592e7x₁⁹x₂¹⁰ + 4.232592e7x₁¹¹x₂⁸ + 2.821728e7x₁¹³x₂⁶ + 1.209312e7x₁¹⁵x₂⁴ + 3.02328e6x₁¹⁷x₂² + 335920.0x₁¹⁹ + 352722.46646x₂²⁰ + 3.5272246646e6x₁²x₂¹⁸ + 1.58725109907e7x₁⁴x₂¹⁶ + 4.23266959752e7x₁⁶x₂¹⁴ + 7.40717179566e7x₁⁸x₂¹² + 8.888606154792e7x₁¹⁰x₂¹⁰ + 7.40717179566e7x₁¹²x₂⁸ + 4.23266959752e7x₁¹⁴x₂⁶ + 1.58725109907e7x₁¹⁶x₂⁴ + 3.5272246646e6x₁¹⁸x₂² + 352722.46646x₁²⁰ + 369512.0x₁x₂²⁰ + 3.69512e6x₁³x₂¹⁸ + 1.662804e7x₁⁵x₂¹⁶ + 4.434144e7x₁⁷x₂¹⁴ + 7.759752e7x₁⁹x₂¹² + 9.3117024e7x₁¹¹x₂¹⁰ + 7.759752e7x₁¹³x₂⁸ + 4.434144e7x₁¹⁵x₂⁶ + 1.662804e7x₁¹⁷x₂⁴ + 3.69512e6x₁¹⁹x₂² + 369512.0x₁²¹ + 352723.05432x₂²² + 3.87995359752e6x₁²x₂²⁰ + 1.93997679876e7x₁⁴x₂¹⁸ + 5.81993039628e7x₁⁶x₂¹⁶ + 1.163986079256e8x₁⁸x₂¹⁴ + 1.6295805109584e8x₁¹⁰x₂¹² + 1.6295805109584e8x₁¹²x₂¹⁰ + 1.163986079256e8x₁¹⁴x₂⁸ + 5.81993039628e7x₁¹⁶x₂⁶ + 1.93997679876e7x₁¹⁸x₂⁴ + 3.87995359752e6x₁²⁰x₂² + 352723.05432x₁²² + 335920.0x₁x₂²² + 3.69512e6x₁³x₂²⁰ + 1.84756e7x₁⁵x₂¹⁸ + 5.54268e7x₁⁷x₂¹⁶ + 1.108536e8x₁⁹x₂¹⁴ + 1.5519504e8x₁¹¹x₂¹² + 1.5519504e8x₁¹³x₂¹⁰ + 1.108536e8x₁¹⁵x₂⁸ + 5.54268e7x₁¹⁷x₂⁶ + 1.84756e7x₁¹⁹x₂⁴ + 3.69512e6x₁²¹x₂² + 335920.0x₁²³ + 293936.46646x₂²⁴ + 3.52723759752e6x₁²x₂²² + 1.939980678636e7x₁⁴x₂²⁰ + 6.46660226212e7x₁⁶x₂¹⁸ + 1.454985508977e8x₁⁸x₂¹⁶ + 2.3279768143632e8x₁¹⁰x₂¹⁴ + 2.7159729500904e8x₁¹²x₂¹² + 2.3279768143632e8x₁¹⁴x₂¹⁰ + 1.454985508977e8x₁¹⁶x₂⁸ + 6.46660226212e7x₁¹⁸x₂⁶ + 1.939980678636e7x₁²⁰x₂⁴ + 3.52723759752e6x₁²²x₂² + 293936.46646x₁²⁴ + 251940.0x₁x₂²⁴ + 3.02328e6x₁³x₂²² + 1.662804e7x₁⁵x₂²⁰ + 5.54268e7x₁⁷x₂¹⁸ + 1.247103e8x₁⁹x₂¹⁶ + 1.9953648e8x₁¹¹x₂¹⁴ + 2.3279256e8x₁¹³x₂¹² + 1.9953648e8x₁¹⁵x₂¹⁰ + 1.247103e8x₁¹⁷x₂⁸ + 5.54268e7x₁¹⁹x₂⁶ + 1.662804e7x₁²¹x₂⁴ + 3.02328e6x₁²³x₂² + 251940.0x₁²⁵ + 203494.9742x₂²⁶ + 2.6454346646e6x₁²x₂²⁴ + 1.58726079876e7x₁⁴x₂²² + 5.81995626212e7x₁⁶x₂²⁰ + 1.45498906553e8x₁⁸x₂¹⁸ + 2.618980317954e8x₁¹⁰x₂¹⁶ + 3.491973757272e8x₁¹²x₂¹⁴ + 3.491973757272e8x₁¹⁴x₂¹² + 2.618980317954e8x₁¹⁶x₂¹⁰ + 1.45498906553e8x₁¹⁸x₂⁸ + 5.81995626212e7x₁²⁰x₂⁶ + 1.58726079876e7x₁²²x₂⁴ + 2.6454346646e6x₁²⁴x₂² + 203494.9742x₁²⁶ + 155040.0x₁x₂²⁶ + 2.01552e6x₁³x₂²⁴ + 1.209312e7x₁⁵x₂²² + 4.434144e7x₁⁷x₂²⁰ + 1.108536e8x₁⁹x₂¹⁸ + 1.9953648e8x₁¹¹x₂¹⁶ + 2.6604864e8x₁¹³x₂¹⁴ + 2.6604864e8x₁¹⁵x₂¹² + 1.9953648e8x₁¹⁷x₂¹⁰ + 1.108536e8x₁¹⁹x₂⁸ + 4.434144e7x₁²¹x₂⁶ + 1.209312e7x₁²³x₂⁴ + 2.01552e6x₁²⁵x₂² + 155040.0x₁²⁷ + 116283.1977x₂²⁸ + 1.6279647678e6x₁²x₂²⁶ + 1.05817709907e7x₁⁴x₂²⁴ + 4.23270839628e7x₁⁶x₂²² + 1.163994808977e8x₁⁸x₂²⁰ + 2.327989617954e8x₁¹⁰x₂¹⁸ + 3.491984426931e8x₁¹²x₂¹⁶ + 3.990839345064e8x₁¹⁴x₂¹⁴ + 3.491984426931e8x₁¹⁶x₂¹² + 2.327989617954e8x₁¹⁸x₂¹⁰ + 1.163994808977e8x₁²⁰x₂⁸ + 4.23270839628e7x₁²²x₂⁶ + 1.05817709907e7x₁²⁴x₂⁴ + 1.6279647678e6x₁²⁶x₂² + 116283.1977x₁²⁸ + 77520.0x₁x₂²⁸ + 1.08528e6x₁³x₂²⁶ + 7.05432e6x₁⁵x₂²⁴ + 2.821728e7x₁⁷x₂²² + 7.759752e7x₁⁹x₂²⁰ + 1.5519504e8x₁¹¹x₂¹⁸ + 2.3279256e8x₁¹³x₂¹⁶ + 2.6604864e8x₁¹⁵x₂¹⁴ + 2.3279256e8x₁¹⁷x₂¹² + 1.5519504e8x₁¹⁹x₂¹⁰ + 7.759752e7x₁²¹x₂⁸ + 2.821728e7x₁²³x₂⁶ + 7.05432e6x₁²⁵x₂⁴ + 1.08528e6x₁²⁷x₂² + 77520.0x₁²⁹ + 54265.70544x₂³⁰ + 813985.5816x₁²x₂²⁸ + 5.6978990712e6x₁⁴x₂²⁶ + 2.46908959752e7x₁⁶x₂²⁴ + 7.40726879256e7x₁⁸x₂²² + 1.6295991343632e8x₁¹⁰x₂²⁰ + 2.715998557272e8x₁¹²x₂¹⁸ + 3.491998145064e8x₁¹⁴x₂¹⁶ + 3.491998145064e8x₁¹⁶x₂¹⁴ + 2.715998557272e8x₁¹⁸x₂¹² + 1.6295991343632e8x₁²⁰x₂¹⁰ + 7.40726879256e7x₁²²x₂⁸ + 2.46908959752e7x₁²⁴x₂⁶ + 5.6978990712e6x₁²⁶x₂⁴ + 813985.5816x₁²⁸x₂² + 54265.70544x₁³⁰ + 31008.0x₁x₂³⁰ + 465120.0x₁³x₂²⁸ + 3.25584e6x₁⁵x₂²⁶ + 1.410864e7x₁⁷x₂²⁴ + 4.232592e7x₁⁹x₂²² + 9.3117024e7x₁¹¹x₂²⁰ + 1.5519504e8x₁¹³x₂¹⁸ + 1.9953648e8x₁¹⁵x₂¹⁶ + 1.9953648e8x₁¹⁷x₂¹⁴ + 1.5519504e8x₁¹⁹x₂¹² + 9.3117024e7x₁²¹x₂¹⁰ + 4.232592e7x₁²³x₂⁸ + 1.410864e7x₁²⁵x₂⁶ + 3.25584e6x₁²⁷x₂⁴ + 465120.0x₁²⁹x₂² + 31008.0x₁³¹ + 20349.74613x₂³² + 325595.93808x₁²x₂³⁰ + 2.4419695356e6x₁⁴x₂²⁸ + 1.13958578328e7x₁⁶x₂²⁶ + 3.70365379566e7x₁⁸x₂²⁴ + 8.888769109584e7x₁¹⁰x₂²² + 1.6296076700904e8x₁¹²x₂²⁰ + 2.328010957272e8x₁¹⁴x₂¹⁸ + 2.619012326931e8x₁¹⁶x₂¹⁶ + 2.328010957272e8x₁¹⁸x₂¹⁴ + 1.6296076700904e8x₁²⁰x₂¹² + 8.888769109584e7x₁²²x₂¹⁰ + 3.70365379566e7x₁²⁴x₂⁸ + 1.13958578328e7x₁²⁶x₂⁶ + 2.4419695356e6x₁²⁸x₂⁴ + 325595.93808x₁³⁰x₂² + 20349.74613x₁³² + 9690.0x₁x₂³² + 155040.0x₁³x₂³⁰ + 1.1628e6x₁⁵x₂²⁸ + 5.4264e6x₁⁷x₂²⁶ + 1.76358e7x₁⁹x₂²⁴ + 4.232592e7x₁¹¹x₂²² + 7.759752e7x₁¹³x₂²⁰ + 1.108536e8x₁¹⁵x₂¹⁸ + 1.247103e8x₁¹⁷x₂¹⁶ + 1.108536e8x₁¹⁹x₂¹⁴ + 7.759752e7x₁²¹x₂¹² + 4.232592e7x₁²³x₂¹⁰ + 1.76358e7x₁²⁵x₂⁸ + 5.4264e6x₁²⁷x₂⁶ + 1.1628e6x₁²⁹x₂⁴ + 155040.0x₁³¹x₂² + 9690.0x₁³³ + 5985.26334x₂³⁴ + 101749.47678x₁²x₂³² + 813995.81424x₁⁴x₂³⁰ + 4.0699790712e6x₁⁶x₂²⁸ + 1.42449267492e7x₁⁸x₂²⁶ + 3.703680954792e7x₁¹⁰x₂²⁴ + 7.407361909584e7x₁¹²x₂²² + 1.1640140143632e8x₁¹⁴x₂²⁰ + 1.455017517954e8x₁¹⁶x₂¹⁸ + 1.455017517954e8x₁¹⁸x₂¹⁶ + 1.1640140143632e8x₁²⁰x₂¹⁴ + 7.407361909584e7x₁²²x₂¹² + 3.703680954792e7x₁²⁴x₂¹⁰ + 1.42449267492e7x₁²⁶x₂⁸ + 4.0699790712e6x₁²⁸x₂⁶ + 813995.81424x₁³⁰x₂⁴ + 101749.47678x₁³²x₂² + 5985.26334x₁³⁴ + 2280.0x₁x₂³⁴ + 38760.0x₁³x₂³² + 310080.0x₁⁵x₂³⁰ + 1.5504e6x₁⁷x₂²⁸ + 5.4264e6x₁⁹x₂²⁶ + 1.410864e7x₁¹¹x₂²⁴ + 2.821728e7x₁¹³x₂²² + 4.434144e7x₁¹⁵x₂²⁰ + 5.54268e7x₁¹⁷x₂¹⁸ + 5.54268e7x₁¹⁹x₂¹⁶ + 4.434144e7x₁²¹x₂¹⁴ + 2.821728e7x₁²³x₂¹² + 1.410864e7x₁²⁵x₂¹⁰ + 5.4264e6x₁²⁷x₂⁸ + 1.5504e6x₁²⁹x₂⁶ + 310080.0x₁³¹x₂⁴ + 38760.0x₁³³x₂² + 2280.0x₁³⁵ + 1330.07315x₂³⁶ + 23941.3167x₁²x₂³⁴ + 203501.19195x₁⁴x₂³² + 1.0853396904e6x₁⁶x₂³⁰ + 4.070023839e6x₁⁸x₂²⁸ + 1.13960667492e7x₁¹⁰x₂²⁶ + 2.46914779566e7x₁¹²x₂²⁴ + 4.23282479256e7x₁¹⁴x₂²² + 5.82013408977e7x₁¹⁶x₂²⁰ + 6.4668156553e7x₁¹⁸x₂¹⁸ + 5.82013408977e7x₁²⁰x₂¹⁶ + 4.23282479256e7x₁²²x₂¹⁴ + 2.46914779566e7x₁²⁴x₂¹² + 1.13960667492e7x₁²⁶x₂¹⁰ + 4.070023839e6x₁²⁸x₂⁸ + 1.0853396904e6x₁³⁰x₂⁶ + 203501.19195x₁³²x₂⁴ + 23941.3167x₁³⁴x₂² + 1330.07315x₁³⁶ + 380.0x₁x₂³⁶ + 6840.0x₁³x₂³⁴ + 58140.0x₁⁵x₂³² + 310080.0x₁⁷x₂³⁰ + 1.1628e6x₁⁹x₂²⁸ + 3.25584e6x₁¹¹x₂²⁶ + 7.05432e6x₁¹³x₂²⁴ + 1.209312e7x₁¹⁵x₂²² + 1.662804e7x₁¹⁷x₂²⁰ + 1.84756e7x₁¹⁹x₂¹⁸ + 1.662804e7x₁²¹x₂¹⁶ + 1.209312e7x₁²³x₂¹⁴ + 7.05432e6x₁²⁵x₂¹² + 3.25584e6x₁²⁷x₂¹⁰ + 1.1628e6x₁²⁹x₂⁸ + 310080.0x₁³¹x₂⁶ + 58140.0x₁³³x₂⁴ + 6840.0x₁³⁵x₂² + 380.0x₁³⁷ + 210.0154x₂³⁸ + 3990.2926x₁²x₂³⁶ + 35912.6334x₁⁴x₂³⁴ + 203504.9226x₁⁶x₂³² + 814019.6904x₁⁸x₂³⁰ + 2.4420590712e6x₁¹⁰x₂²⁸ + 5.6981378328e6x₁¹²x₂²⁶ + 1.05822559752e7x₁¹⁴x₂²⁴ + 1.58733839628e7x₁¹⁶x₂²² + 1.94008026212e7x₁¹⁸x₂²⁰ + 1.94008026212e7x₁²⁰x₂¹⁸ + 1.58733839628e7x₁²²x₂¹⁶ + 1.05822559752e7x₁²⁴x₂¹⁴ + 5.6981378328e6x₁²⁶x₂¹² + 2.4420590712e6x₁²⁸x₂¹⁰ + 814019.6904x₁³⁰x₂⁸ + 203504.9226x₁³²x₂⁶ + 35912.6334x₁³⁴x₂⁴ + 3990.2926x₁³⁶x₂² + 210.0154x₁³⁸ + 40.0x₁x₂³⁸ + 760.0x₁³x₂³⁶ + 6840.0x₁⁵x₂³⁴ + 38760.0x₁⁷x₂³² + 155040.0x₁⁹x₂³⁰ + 465120.0x₁¹¹x₂²⁸ + 1.08528e6x₁¹³x₂²⁶ + 2.01552e6x₁¹⁵x₂²⁴ + 3.02328e6x₁¹⁷x₂²² + 3.69512e6x₁¹⁹x₂²⁰ + 3.69512e6x₁²¹x₂¹⁸ + 3.02328e6x₁²³x₂¹⁶ + 2.01552e6x₁²⁵x₂¹⁴ + 1.08528e6x₁²⁷x₂¹² + 465120.0x₁²⁹x₂¹⁰ + 155040.0x₁³¹x₂⁸ + 38760.0x₁³³x₂⁶ + 6840.0x₁³⁵x₂⁴ + 760.0x₁³⁷x₂² + 40.0x₁³⁹ + 21.00231x₂⁴⁰ + 420.0462x₁²x₂³⁸ + 3990.4389x₁⁴x₂³⁶ + 23942.6334x₁⁶x₂³⁴ + 101756.19195x₁⁸x₂³² + 325619.81424x₁¹⁰x₂³⁰ + 814049.5356x₁¹²x₂²⁸ + 1.6280990712e6x₁¹⁴x₂²⁶ + 2.6456609907e6x₁¹⁶x₂²⁴ + 3.5275479876e6x₁¹⁸x₂²² + 3.88030278636e6x₁²⁰x₂²⁰ + 3.5275479876e6x₁²²x₂¹⁸ + 2.6456609907e6x₁²⁴x₂¹⁶ + 1.6280990712e6x₁²⁶x₂¹⁴ + 814049.5356x₁²⁸x₂¹² + 325619.81424x₁³⁰x₂¹⁰ + 101756.19195x₁³²x₂⁸ + 23942.6334x₁³⁴x₂⁶ + 3990.4389x₁³⁶x₂⁴ + 420.0462x₁³⁸x₂² + 21.00231x₁⁴⁰ + 2.0x₁x₂⁴⁰ + 40.0x₁³x₂³⁸ + 380.0x₁⁵x₂³⁶ + 2280.0x₁⁷x₂³⁴ + 9690.0x₁⁹x₂³² + 31008.0x₁¹¹x₂³⁰ + 77520.0x₁¹³x₂²⁸ + 155040.0x₁¹⁵x₂²⁶ + 251940.0x₁¹⁷x₂²⁴ + 335920.0x₁¹⁹x₂²² + 369512.0x₁²¹x₂²⁰ + 335920.0x₁²³x₂¹⁸ + 251940.0x₁²⁵x₂¹⁶ + 155040.0x₁²⁷x₂¹⁴ + 77520.0x₁²⁹x₂¹² + 31008.0x₁³¹x₂¹⁰ + 9690.0x₁³³x₂⁸ + 2280.0x₁³⁵x₂⁶ + 380.0x₁³⁷x₂⁴ + 40.0x₁³⁹x₂² + 2.0x₁⁴¹ + 1.00022x₂⁴² + 21.00462x₁²x₂⁴⁰ + 210.0462x₁⁴x₂³⁸ + 1330.2926x₁⁶x₂³⁶ + 5986.3167x₁⁸x₂³⁴ + 20353.47678x₁¹⁰x₂³² + 54275.93808x₁¹²x₂³⁰ + 116305.5816x₁¹⁴x₂²⁸ + 203534.7678x₁¹⁶x₂²⁶ + 293994.6646x₁¹⁸x₂²⁴ + 352793.59752x₁²⁰x₂²² + 352793.59752x₁²²x₂²⁰ + 293994.6646x₁²⁴x₂¹⁸ + 203534.7678x₁²⁶x₂¹⁶ + 116305.5816x₁²⁸x₂¹⁴ + 54275.93808x₁³⁰x₂¹² + 20353.47678x₁³²x₂¹⁰ + 5986.3167x₁³⁴x₂⁸ + 1330.2926x₁³⁶x₂⁶ + 210.0462x₁³⁸x₂⁴ + 21.00462x₁⁴⁰x₂² + 1.00022x₁⁴² + 1.0e-5x₂⁴⁴ + 0.00022x₁²x₂⁴² + 0.00231x₁⁴x₂⁴⁰ + 0.0154x₁⁶x₂³⁸ + 0.07315x₁⁸x₂³⁶ + 0.26334x₁¹⁰x₂³⁴ + 0.74613x₁¹²x₂³² + 1.70544x₁¹⁴x₂³⁰ + 3.1977x₁¹⁶x₂²⁸ + 4.9742x₁¹⁸x₂²⁶ + 6.46646x₁²⁰x₂²⁴ + 7.05432x₁²²x₂²² + 6.46646x₁²⁴x₂²⁰ + 4.9742x₁²⁶x₂¹⁸ + 3.1977x₁²⁸x₂¹⁶ + 1.70544x₁³⁰x₂¹⁴ + 0.74613x₁³²x₂¹² + 0.26334x₁³⁴x₂¹⁰ + 0.07315x₁³⁶x₂⁸ + 0.0154x₁³⁸x₂⁶ + 0.00231x₁⁴⁰x₂⁴ + 0.00022x₁⁴²x₂² + 1.0e-5x₁⁴⁴
Objective was scaled by the prefactor 1.0 + 20.0x₂² + 20.0x₁² + 190.0x₂⁴ + 380.0x₁²x₂² + 190.0x₁⁴ + 1140.0x₂⁶ + 3420.0x₁²x₂⁴ + 3420.0x₁⁴x₂² + 1140.0x₁⁶ + 4845.0x₂⁸ + 19380.0x₁²x₂⁶ + 29070.0x₁⁴x₂⁴ + 19380.0x₁⁶x₂² + 4845.0x₁⁸ + 15504.0x₂¹⁰ + 77520.0x₁²x₂⁸ + 155040.0x₁⁴x₂⁶ + 155040.0x₁⁶x₂⁴ + 77520.0x₁⁸x₂² + 15504.0x₁¹⁰ + 38760.0x₂¹² + 232560.0x₁²x₂¹⁰ + 581400.0x₁⁴x₂⁸ + 775200.0x₁⁶x₂⁶ + 581400.0x₁⁸x₂⁴ + 232560.0x₁¹⁰x₂² + 38760.0x₁¹² + 77520.0x₂¹⁴ + 542640.0x₁²x₂¹² + 1.62792e6x₁⁴x₂¹⁰ + 2.7132e6x₁⁶x₂⁸ + 2.7132e6x₁⁸x₂⁶ + 1.62792e6x₁¹⁰x₂⁴ + 542640.0x₁¹²x₂² + 77520.0x₁¹⁴ + 125970.0x₂¹⁶ + 1.00776e6x₁²x₂¹⁴ + 3.52716e6x₁⁴x₂¹² + 7.05432e6x₁⁶x₂¹⁰ + 8.8179e6x₁⁸x₂⁸ + 7.05432e6x₁¹⁰x₂⁶ + 3.52716e6x₁¹²x₂⁴ + 1.00776e6x₁¹⁴x₂² + 125970.0x₁¹⁶ + 167960.0x₂¹⁸ + 1.51164e6x₁²x₂¹⁶ + 6.04656e6x₁⁴x₂¹⁴ + 1.410864e7x₁⁶x₂¹² + 2.116296e7x₁⁸x₂¹⁰ + 2.116296e7x₁¹⁰x₂⁸ + 1.410864e7x₁¹²x₂⁶ + 6.04656e6x₁¹⁴x₂⁴ + 1.51164e6x₁¹⁶x₂² + 167960.0x₁¹⁸ + 184756.0x₂²⁰ + 1.84756e6x₁²x₂¹⁸ + 8.31402e6x₁⁴x₂¹⁶ + 2.217072e7x₁⁶x₂¹⁴ + 3.879876e7x₁⁸x₂¹² + 4.6558512e7x₁¹⁰x₂¹⁰ + 3.879876e7x₁¹²x₂⁸ + 2.217072e7x₁¹⁴x₂⁶ + 8.31402e6x₁¹⁶x₂⁴ + 1.84756e6x₁¹⁸x₂² + 184756.0x₁²⁰ + 167960.0x₂²² + 1.84756e6x₁²x₂²⁰ + 9.2378e6x₁⁴x₂¹⁸ + 2.77134e7x₁⁶x₂¹⁶ + 5.54268e7x₁⁸x₂¹⁴ + 7.759752e7x₁¹⁰x₂¹² + 7.759752e7x₁¹²x₂¹⁰ + 5.54268e7x₁¹⁴x₂⁸ + 2.77134e7x₁¹⁶x₂⁶ + 9.2378e6x₁¹⁸x₂⁴ + 1.84756e6x₁²⁰x₂² + 167960.0x₁²² + 125970.0x₂²⁴ + 1.51164e6x₁²x₂²² + 8.31402e6x₁⁴x₂²⁰ + 2.77134e7x₁⁶x₂¹⁸ + 6.235515e7x₁⁸x₂¹⁶ + 9.976824e7x₁¹⁰x₂¹⁴ + 1.1639628e8x₁¹²x₂¹² + 9.976824e7x₁¹⁴x₂¹⁰ + 6.235515e7x₁¹⁶x₂⁸ + 2.77134e7x₁¹⁸x₂⁶ + 8.31402e6x₁²⁰x₂⁴ + 1.51164e6x₁²²x₂² + 125970.0x₁²⁴ + 77520.0x₂²⁶ + 1.00776e6x₁²x₂²⁴ + 6.04656e6x₁⁴x₂²² + 2.217072e7x₁⁶x₂²⁰ + 5.54268e7x₁⁸x₂¹⁸ + 9.976824e7x₁¹⁰x₂¹⁶ + 1.3302432e8x₁¹²x₂¹⁴ + 1.3302432e8x₁¹⁴x₂¹² + 9.976824e7x₁¹⁶x₂¹⁰ + 5.54268e7x₁¹⁸x₂⁸ + 2.217072e7x₁²⁰x₂⁶ + 6.04656e6x₁²²x₂⁴ + 1.00776e6x₁²⁴x₂² + 77520.0x₁²⁶ + 38760.0x₂²⁸ + 542640.0x₁²x₂²⁶ + 3.52716e6x₁⁴x₂²⁴ + 1.410864e7x₁⁶x₂²² + 3.879876e7x₁⁸x₂²⁰ + 7.759752e7x₁¹⁰x₂¹⁸ + 1.1639628e8x₁¹²x₂¹⁶ + 1.3302432e8x₁¹⁴x₂¹⁴ + 1.1639628e8x₁¹⁶x₂¹² + 7.759752e7x₁¹⁸x₂¹⁰ + 3.879876e7x₁²⁰x₂⁸ + 1.410864e7x₁²²x₂⁶ + 3.52716e6x₁²⁴x₂⁴ + 542640.0x₁²⁶x₂² + 38760.0x₁²⁸ + 15504.0x₂³⁰ + 232560.0x₁²x₂²⁸ + 1.62792e6x₁⁴x₂²⁶ + 7.05432e6x₁⁶x₂²⁴ + 2.116296e7x₁⁸x₂²² + 4.6558512e7x₁¹⁰x₂²⁰ + 7.759752e7x₁¹²x₂¹⁸ + 9.976824e7x₁¹⁴x₂¹⁶ + 9.976824e7x₁¹⁶x₂¹⁴ + 7.759752e7x₁¹⁸x₂¹² + 4.6558512e7x₁²⁰x₂¹⁰ + 2.116296e7x₁²²x₂⁸ + 7.05432e6x₁²⁴x₂⁶ + 1.62792e6x₁²⁶x₂⁴ + 232560.0x₁²⁸x₂² + 15504.0x₁³⁰ + 4845.0x₂³² + 77520.0x₁²x₂³⁰ + 581400.0x₁⁴x₂²⁸ + 2.7132e6x₁⁶x₂²⁶ + 8.8179e6x₁⁸x₂²⁴ + 2.116296e7x₁¹⁰x₂²² + 3.879876e7x₁¹²x₂²⁰ + 5.54268e7x₁¹⁴x₂¹⁸ + 6.235515e7x₁¹⁶x₂¹⁶ + 5.54268e7x₁¹⁸x₂¹⁴ + 3.879876e7x₁²⁰x₂¹² + 2.116296e7x₁²²x₂¹⁰ + 8.8179e6x₁²⁴x₂⁸ + 2.7132e6x₁²⁶x₂⁶ + 581400.0x₁²⁸x₂⁴ + 77520.0x₁³⁰x₂² + 4845.0x₁³² + 1140.0x₂³⁴ + 19380.0x₁²x₂³² + 155040.0x₁⁴x₂³⁰ + 775200.0x₁⁶x₂²⁸ + 2.7132e6x₁⁸x₂²⁶ + 7.05432e6x₁¹⁰x₂²⁴ + 1.410864e7x₁¹²x₂²² + 2.217072e7x₁¹⁴x₂²⁰ + 2.77134e7x₁¹⁶x₂¹⁸ + 2.77134e7x₁¹⁸x₂¹⁶ + 2.217072e7x₁²⁰x₂¹⁴ + 1.410864e7x₁²²x₂¹² + 7.05432e6x₁²⁴x₂¹⁰ + 2.7132e6x₁²⁶x₂⁸ + 775200.0x₁²⁸x₂⁶ + 155040.0x₁³⁰x₂⁴ + 19380.0x₁³²x₂² + 1140.0x₁³⁴ + 190.0x₂³⁶ + 3420.0x₁²x₂³⁴ + 29070.0x₁⁴x₂³² + 155040.0x₁⁶x₂³⁰ + 581400.0x₁⁸x₂²⁸ + 1.62792e6x₁¹⁰x₂²⁶ + 3.52716e6x₁¹²x₂²⁴ + 6.04656e6x₁¹⁴x₂²² + 8.31402e6x₁¹⁶x₂²⁰ + 9.2378e6x₁¹⁸x₂¹⁸ + 8.31402e6x₁²⁰x₂¹⁶ + 6.04656e6x₁²²x₂¹⁴ + 3.52716e6x₁²⁴x₂¹² + 1.62792e6x₁²⁶x₂¹⁰ + 581400.0x₁²⁸x₂⁸ + 155040.0x₁³⁰x₂⁶ + 29070.0x₁³²x₂⁴ + 3420.0x₁³⁴x₂² + 190.0x₁³⁶ + 20.0x₂³⁸ + 380.0x₁²x₂³⁶ + 3420.0x₁⁴x₂³⁴ + 19380.0x₁⁶x₂³² + 77520.0x₁⁸x₂³⁰ + 232560.0x₁¹⁰x₂²⁸ + 542640.0x₁¹²x₂²⁶ + 1.00776e6x₁¹⁴x₂²⁴ + 1.51164e6x₁¹⁶x₂²² + 1.84756e6x₁¹⁸x₂²⁰ + 1.84756e6x₁²⁰x₂¹⁸ + 1.51164e6x₁²²x₂¹⁶ + 1.00776e6x₁²⁴x₂¹⁴ + 542640.0x₁²⁶x₂¹² + 232560.0x₁²⁸x₂¹⁰ + 77520.0x₁³⁰x₂⁸ + 19380.0x₁³²x₂⁶ + 3420.0x₁³⁴x₂⁴ + 380.0x₁³⁶x₂² + 20.0x₁³⁸ + x₂⁴⁰ + 20.0x₁²x₂³⁸ + 190.0x₁⁴x₂³⁶ + 1140.0x₁⁶x₂³⁴ + 4845.0x₁⁸x₂³² + 15504.0x₁¹⁰x₂³⁰ + 38760.0x₁¹²x₂²⁸ + 77520.0x₁¹⁴x₂²⁶ + 125970.0x₁¹⁶x₂²⁴ + 167960.0x₁¹⁸x₂²² + 184756.0x₁²⁰x₂²⁰ + 167960.0x₁²²x₂¹⁸ + 125970.0x₁²⁴x₂¹⁶ + 77520.0x₁²⁶x₂¹⁴ + 38760.0x₁²⁸x₂¹² + 15504.0x₁³⁰x₂¹⁰ + 4845.0x₁³²x₂⁸ + 1140.0x₁³⁴x₂⁶ + 190.0x₁³⁶x₂⁴ + 20.0x₁³⁸x₂² + x₁⁴⁰
1 equality constraint
1: -x₂² + x₁³ = 0"
    # this cannot be solved by any. The problem spans too many orders of magnitude anyway.
end

@testset "POP 17 (Example A.8)" begin
    DynamicPolynomials.@polyvar x[1:10]
    Random.seed!(1234)
    prob = poly_problem(sum(x[j]^2 + x[j+5]^2 for j in 1:5)/6,
        zero=[x[6]-1; [x[j+6]-x[j+5]-(x[j+5]^2-x[j])/6 for j in 1:4]], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 10 variables
Objective: 1.0e-5 + 0.1667x₁₀² + 0.1667x₉² + 0.1667x₈² + 0.1667x₇² + 0.1667x₆² + 0.1667x₅² + 0.1667x₄² + 0.1667x₃² + 0.1667x₂² + 0.1667x₁² + 0.1667x₁₀⁴ + 0.33339x₉²x₁₀² + 0.1667x₉⁴ + 0.33339x₈²x₁₀² + 0.33339x₈²x₉² + 0.1667x₈⁴ + 0.33339x₇²x₁₀² + 0.33339x₇²x₉² + 0.33339x₇²x₈² + 0.1667x₇⁴ + 0.33339x₆²x₁₀² + 0.33339x₆²x₉² + 0.33339x₆²x₈² + 0.33339x₆²x₇² + 0.1667x₆⁴ + 0.33339x₅²x₁₀² + 0.33339x₅²x₉² + 0.33339x₅²x₈² + 0.33339x₅²x₇² + 0.33339x₅²x₆² + 0.1667x₅⁴ + 0.33339x₄²x₁₀² + 0.33339x₄²x₉² + 0.33339x₄²x₈² + 0.33339x₄²x₇² + 0.33339x₄²x₆² + 0.33339x₄²x₅² + 0.1667x₄⁴ + 0.33339x₃²x₁₀² + 0.33339x₃²x₉² + 0.33339x₃²x₈² + 0.33339x₃²x₇² + 0.33339x₃²x₆² + 0.33339x₃²x₅² + 0.33339x₃²x₄² + 0.1667x₃⁴ + 0.33339x₂²x₁₀² + 0.33339x₂²x₉² + 0.33339x₂²x₈² + 0.33339x₂²x₇² + 0.33339x₂²x₆² + 0.33339x₂²x₅² + 0.33339x₂²x₄² + 0.33339x₂²x₃² + 0.1667x₂⁴ + 0.33339x₁²x₁₀² + 0.33339x₁²x₉² + 0.33339x₁²x₈² + 0.33339x₁²x₇² + 0.33339x₁²x₆² + 0.33339x₁²x₅² + 0.33339x₁²x₄² + 0.33339x₁²x₃² + 0.33339x₁²x₂² + 0.1667x₁⁴ + 1.0e-5x₁₀⁶ + 3.0e-5x₉²x₁₀⁴ + 3.0e-5x₉⁴x₁₀² + 1.0e-5x₉⁶ + 3.0e-5x₈²x₁₀⁴ + 6.0e-5x₈²x₉²x₁₀² + 3.0e-5x₈²x₉⁴ + 3.0e-5x₈⁴x₁₀² + 3.0e-5x₈⁴x₉² + 1.0e-5x₈⁶ + 3.0e-5x₇²x₁₀⁴ + 6.0e-5x₇²x₉²x₁₀² + 3.0e-5x₇²x₉⁴ + 6.0e-5x₇²x₈²x₁₀² + 6.0e-5x₇²x₈²x₉² + 3.0e-5x₇²x₈⁴ + 3.0e-5x₇⁴x₁₀² + 3.0e-5x₇⁴x₉² + 3.0e-5x₇⁴x₈² + 1.0e-5x₇⁶ + 3.0e-5x₆²x₁₀⁴ + 6.0e-5x₆²x₉²x₁₀² + 3.0e-5x₆²x₉⁴ + 6.0e-5x₆²x₈²x₁₀² + 6.0e-5x₆²x₈²x₉² + 3.0e-5x₆²x₈⁴ + 6.0e-5x₆²x₇²x₁₀² + 6.0e-5x₆²x₇²x₉² + 6.0e-5x₆²x₇²x₈² + 3.0e-5x₆²x₇⁴ + 3.0e-5x₆⁴x₁₀² + 3.0e-5x₆⁴x₉² + 3.0e-5x₆⁴x₈² + 3.0e-5x₆⁴x₇² + 1.0e-5x₆⁶ + 3.0e-5x₅²x₁₀⁴ + 6.0e-5x₅²x₉²x₁₀² + 3.0e-5x₅²x₉⁴ + 6.0e-5x₅²x₈²x₁₀² + 6.0e-5x₅²x₈²x₉² + 3.0e-5x₅²x₈⁴ + 6.0e-5x₅²x₇²x₁₀² + 6.0e-5x₅²x₇²x₉² + 6.0e-5x₅²x₇²x₈² + 3.0e-5x₅²x₇⁴ + 6.0e-5x₅²x₆²x₁₀² + 6.0e-5x₅²x₆²x₉² + 6.0e-5x₅²x₆²x₈² + 6.0e-5x₅²x₆²x₇² + 3.0e-5x₅²x₆⁴ + 3.0e-5x₅⁴x₁₀² + 3.0e-5x₅⁴x₉² + 3.0e-5x₅⁴x₈² + 3.0e-5x₅⁴x₇² + 3.0e-5x₅⁴x₆² + 1.0e-5x₅⁶ + 3.0e-5x₄²x₁₀⁴ + 6.0e-5x₄²x₉²x₁₀² + 3.0e-5x₄²x₉⁴ + 6.0e-5x₄²x₈²x₁₀² + 6.0e-5x₄²x₈²x₉² + 3.0e-5x₄²x₈⁴ + 6.0e-5x₄²x₇²x₁₀² + 6.0e-5x₄²x₇²x₉² + 6.0e-5x₄²x₇²x₈² + 3.0e-5x₄²x₇⁴ + 6.0e-5x₄²x₆²x₁₀² + 6.0e-5x₄²x₆²x₉² + 6.0e-5x₄²x₆²x₈² + 6.0e-5x₄²x₆²x₇² + 3.0e-5x₄²x₆⁴ + 6.0e-5x₄²x₅²x₁₀² + 6.0e-5x₄²x₅²x₉² + 6.0e-5x₄²x₅²x₈² + 6.0e-5x₄²x₅²x₇² + 6.0e-5x₄²x₅²x₆² + 3.0e-5x₄²x₅⁴ + 3.0e-5x₄⁴x₁₀² + 3.0e-5x₄⁴x₉² + 3.0e-5x₄⁴x₈² + 3.0e-5x₄⁴x₇² + 3.0e-5x₄⁴x₆² + 3.0e-5x₄⁴x₅² + 1.0e-5x₄⁶ + 3.0e-5x₃²x₁₀⁴ + 6.0e-5x₃²x₉²x₁₀² + 3.0e-5x₃²x₉⁴ + 6.0e-5x₃²x₈²x₁₀² + 6.0e-5x₃²x₈²x₉² + 3.0e-5x₃²x₈⁴ + 6.0e-5x₃²x₇²x₁₀² + 6.0e-5x₃²x₇²x₉² + 6.0e-5x₃²x₇²x₈² + 3.0e-5x₃²x₇⁴ + 6.0e-5x₃²x₆²x₁₀² + 6.0e-5x₃²x₆²x₉² + 6.0e-5x₃²x₆²x₈² + 6.0e-5x₃²x₆²x₇² + 3.0e-5x₃²x₆⁴ + 6.0e-5x₃²x₅²x₁₀² + 6.0e-5x₃²x₅²x₉² + 6.0e-5x₃²x₅²x₈² + 6.0e-5x₃²x₅²x₇² + 6.0e-5x₃²x₅²x₆² + 3.0e-5x₃²x₅⁴ + 6.0e-5x₃²x₄²x₁₀² + 6.0e-5x₃²x₄²x₉² + 6.0e-5x₃²x₄²x₈² + 6.0e-5x₃²x₄²x₇² + 6.0e-5x₃²x₄²x₆² + 6.0e-5x₃²x₄²x₅² + 3.0e-5x₃²x₄⁴ + 3.0e-5x₃⁴x₁₀² + 3.0e-5x₃⁴x₉² + 3.0e-5x₃⁴x₈² + 3.0e-5x₃⁴x₇² + 3.0e-5x₃⁴x₆² + 3.0e-5x₃⁴x₅² + 3.0e-5x₃⁴x₄² + 1.0e-5x₃⁶ + 3.0e-5x₂²x₁₀⁴ + 6.0e-5x₂²x₉²x₁₀² + 3.0e-5x₂²x₉⁴ + 6.0e-5x₂²x₈²x₁₀² + 6.0e-5x₂²x₈²x₉² + 3.0e-5x₂²x₈⁴ + 6.0e-5x₂²x₇²x₁₀² + 6.0e-5x₂²x₇²x₉² + 6.0e-5x₂²x₇²x₈² + 3.0e-5x₂²x₇⁴ + 6.0e-5x₂²x₆²x₁₀² + 6.0e-5x₂²x₆²x₉² + 6.0e-5x₂²x₆²x₈² + 6.0e-5x₂²x₆²x₇² + 3.0e-5x₂²x₆⁴ + 6.0e-5x₂²x₅²x₁₀² + 6.0e-5x₂²x₅²x₉² + 6.0e-5x₂²x₅²x₈² + 6.0e-5x₂²x₅²x₇² + 6.0e-5x₂²x₅²x₆² + 3.0e-5x₂²x₅⁴ + 6.0e-5x₂²x₄²x₁₀² + 6.0e-5x₂²x₄²x₉² + 6.0e-5x₂²x₄²x₈² + 6.0e-5x₂²x₄²x₇² + 6.0e-5x₂²x₄²x₆² + 6.0e-5x₂²x₄²x₅² + 3.0e-5x₂²x₄⁴ + 6.0e-5x₂²x₃²x₁₀² + 6.0e-5x₂²x₃²x₉² + 6.0e-5x₂²x₃²x₈² + 6.0e-5x₂²x₃²x₇² + 6.0e-5x₂²x₃²x₆² + 6.0e-5x₂²x₃²x₅² + 6.0e-5x₂²x₃²x₄² + 3.0e-5x₂²x₃⁴ + 3.0e-5x₂⁴x₁₀² + 3.0e-5x₂⁴x₉² + 3.0e-5x₂⁴x₈² + 3.0e-5x₂⁴x₇² + 3.0e-5x₂⁴x₆² + 3.0e-5x₂⁴x₅² + 3.0e-5x₂⁴x₄² + 3.0e-5x₂⁴x₃² + 1.0e-5x₂⁶ + 3.0e-5x₁²x₁₀⁴ + 6.0e-5x₁²x₉²x₁₀² + 3.0e-5x₁²x₉⁴ + 6.0e-5x₁²x₈²x₁₀² + 6.0e-5x₁²x₈²x₉² + 3.0e-5x₁²x₈⁴ + 6.0e-5x₁²x₇²x₁₀² + 6.0e-5x₁²x₇²x₉² + 6.0e-5x₁²x₇²x₈² + 3.0e-5x₁²x₇⁴ + 6.0e-5x₁²x₆²x₁₀² + 6.0e-5x₁²x₆²x₉² + 6.0e-5x₁²x₆²x₈² + 6.0e-5x₁²x₆²x₇² + 3.0e-5x₁²x₆⁴ + 6.0e-5x₁²x₅²x₁₀² + 6.0e-5x₁²x₅²x₉² + 6.0e-5x₁²x₅²x₈² + 6.0e-5x₁²x₅²x₇² + 6.0e-5x₁²x₅²x₆² + 3.0e-5x₁²x₅⁴ + 6.0e-5x₁²x₄²x₁₀² + 6.0e-5x₁²x₄²x₉² + 6.0e-5x₁²x₄²x₈² + 6.0e-5x₁²x₄²x₇² + 6.0e-5x₁²x₄²x₆² + 6.0e-5x₁²x₄²x₅² + 3.0e-5x₁²x₄⁴ + 6.0e-5x₁²x₃²x₁₀² + 6.0e-5x₁²x₃²x₉² + 6.0e-5x₁²x₃²x₈² + 6.0e-5x₁²x₃²x₇² + 6.0e-5x₁²x₃²x₆² + 6.0e-5x₁²x₃²x₅² + 6.0e-5x₁²x₃²x₄² + 3.0e-5x₁²x₃⁴ + 6.0e-5x₁²x₂²x₁₀² + 6.0e-5x₁²x₂²x₉² + 6.0e-5x₁²x₂²x₈² + 6.0e-5x₁²x₂²x₇² + 6.0e-5x₁²x₂²x₆² + 6.0e-5x₁²x₂²x₅² + 6.0e-5x₁²x₂²x₄² + 6.0e-5x₁²x₂²x₃² + 3.0e-5x₁²x₂⁴ + 3.0e-5x₁⁴x₁₀² + 3.0e-5x₁⁴x₉² + 3.0e-5x₁⁴x₈² + 3.0e-5x₁⁴x₇² + 3.0e-5x₁⁴x₆² + 3.0e-5x₁⁴x₅² + 3.0e-5x₁⁴x₄² + 3.0e-5x₁⁴x₃² + 3.0e-5x₁⁴x₂² + 1.0e-5x₁⁶
Objective was scaled by the prefactor 1.0 + x₁₀² + x₉² + x₈² + x₇² + x₆² + x₅² + x₄² + x₃² + x₂² + x₁²
5 equality constraints
1: -1.0 + x₆ = 0
2: x₇ - x₆ + 0.16666666666666666x₁ - 0.16666666666666666x₆² = 0
3: x₈ - x₇ + 0.16666666666666666x₂ - 0.16666666666666666x₇² = 0
4: x₉ - x₈ + 0.16666666666666666x₃ - 0.16666666666666666x₈² = 0
5: x₁₀ - x₉ + 0.16666666666666666x₄ - 0.16666666666666666x₉² = 0"
    if optimize
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, prob).objective ≈ 1.321664 atol = 1e-6
        #:COPTMoment ∈ solvers && @test poly_optimize(:COPTMoment, prob).objective ≈ 1.321664 atol = 1e-6 # 70 seconds
    end
end

@testset "POP 18 (Self made)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(x[1]^6 + x[2]^2, zero=[(x[1]^2 + x[2]^2)*(1 - x[1]*x[2])^2], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 + 1.00005x₂² + 5.0e-5x₁² + 1.0001x₂⁴ + 1.0002x₁²x₂² + 0.0001x₁⁴ + 0.0001x₂⁶ + 0.0003x₁²x₂⁴ + 0.0003x₁⁴x₂² + 1.0001x₁⁶ + 5.0e-5x₂⁸ + 0.0002x₁²x₂⁶ + 0.0003x₁⁴x₂⁴ + 1.0002x₁⁶x₂² + 1.00005x₁⁸ + 1.0e-5x₂¹⁰ + 5.0e-5x₁²x₂⁸ + 0.0001x₁⁴x₂⁶ + 0.0001x₁⁶x₂⁴ + 5.0e-5x₁⁸x₂² + 1.0e-5x₁¹⁰
Objective was scaled by the prefactor 1.0 + x₂² + x₁²
1 equality constraint
1: x₂² + x₁² - 2.0x₁x₂³ - 2.0x₁³x₂ + x₁²x₂⁴ + x₁⁴x₂² = 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0 atol = 4e-5
            end
        end
    end
end

@testset "POP 19 (Example 2)" begin
    DynamicPolynomials.@polyvar x[1:4]
    Random.seed!(1234)
    prob = poly_problem(1.0sum(x .^ 2), nonneg=[1/8-x[4]], zero=[sum(x)-1], noncompact=(1e-5, 0))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 4 variables
Objective: 1.0e-5 + 1.00002x₄² + 1.00002x₃² + 1.00002x₂² + 1.00002x₁² + 1.0e-5x₄⁴ + 2.0e-5x₃²x₄² + 1.0e-5x₃⁴ + 2.0e-5x₂²x₄² + 2.0e-5x₂²x₃² + 1.0e-5x₂⁴ + 2.0e-5x₁²x₄² + 2.0e-5x₁²x₃² + 2.0e-5x₁²x₂² + 1.0e-5x₁⁴
1 equality constraint
1: -1.0 + x₄ + x₃ + x₂ + x₁ = 0
1 nonnegative constraint
2: 0.125 - x₄ ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 0.2708494 atol = 2e-4
            end
        end
    end
end

@testset "POP 20 (Self made)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(x[1]^3 - x[2]^2, nonneg=x, zero=[(x[1]*x[2] +1)*(x[1] - x[2])^2], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 - 0.99997x₂² + 3.0e-5x₁² + x₁³ - 0.99997x₂⁴ - 0.99994x₁²x₂² + 3.0e-5x₁⁴ + x₁³x₂² + x₁⁵ + 1.0e-5x₂⁶ + 3.0e-5x₁²x₂⁴ + 3.0e-5x₁⁴x₂² + 1.0e-5x₁⁶
Objective was scaled by the prefactor 1.0 + x₂² + x₁²
1 equality constraint
1: x₂² - 2.0x₁x₂ + x₁² + x₁x₂³ - 2.0x₁²x₂² + x₁³x₂ = 0
2 nonnegative constraints
2: x₁ ≥ 0
3: x₂ ≥ 0"
    # note that this problem is not well-defined. All solvers claim that they don't find the optimal solution (although for
    # the interior-points, it's pretty close).
end

@testset "POP 21 (Self made)" begin
    DynamicPolynomials.@polyvar x[1:2]
    Random.seed!(1234)
    prob = poly_problem(x[1]^4 - 3x[2], nonneg=x, zero=[(x[2] - x[1]^2)*(2x[1]^2 - x[2])], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 2 variables
Objective: 1.0e-5 - 3.0x₂ + 4.0e-5x₂² + 4.0e-5x₁² - 3.0x₂³ - 3.0x₁²x₂ + 6.0e-5x₂⁴ + 0.00012x₁²x₂² + 1.00006x₁⁴ + 4.0e-5x₂⁶ + 0.00012x₁²x₂⁴ + 1.00012x₁⁴x₂² + 1.00004x₁⁶ + 1.0e-5x₂⁸ + 4.0e-5x₁²x₂⁶ + 6.0e-5x₁⁴x₂⁴ + 4.0e-5x₁⁶x₂² + 1.0e-5x₁⁸
Objective was scaled by the prefactor 1.0 + x₂² + x₁²
1 equality constraint
1: -x₂² + 3.0x₁²x₂ - 2.0x₁⁴ = 0
2 nonnegative constraints
2: x₁ ≥ 0
3: x₂ ≥ 0"
    if optimize
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, prob).objective ≈ -8.5578 atol = 1e-4
        :COPTSOS ∈ solvers && @test poly_optimize(:COPTSOS, prob).objective ≈ -8.5578 atol = 1e-4
        :HypatiaMoment ∈ solvers && @test poly_optimize(:HypatiaMoment, prob).objective ≈ -8.5578 atol = 1e-4
        :SCSMoment ∈ solvers && @test poly_optimize(:SCSMoment, prob).objective ≈ -8.5578 atol = 1e-4 broken = true
        # MosekMoment: UNKNOWN, Clarabel: NUMERICAL_ERROR, SCS claims good status, although the result is complete off (-7)
    end
end

@testset "POP 22 (AM-GM inequality)" begin
    DynamicPolynomials.@polyvar x[1:3]
    Random.seed!(1234)
    prob = poly_problem(sum(x), nonneg=x, zero=[prod(x)-1], noncompact=(1e-5, 2))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: 1.0e-5 + x₃ + x₂ + x₁ + 3.0e-5x₃² + 3.0e-5x₂² + 3.0e-5x₁² + 2.0x₃³ + 2.0x₂x₃² + 2.0x₂²x₃ + 2.0x₂³ + 2.0x₁x₃² + 2.0x₁x₂² + 2.0x₁²x₃ + 2.0x₁²x₂ + 2.0x₁³ + 3.0e-5x₃⁴ + 6.0e-5x₂²x₃² + 3.0e-5x₂⁴ + 6.0e-5x₁²x₃² + 6.0e-5x₁²x₂² + 3.0e-5x₁⁴ + x₃⁵ + x₂x₃⁴ + 2.0x₂²x₃³ + 2.0x₂³x₃² + x₂⁴x₃ + x₂⁵ + x₁x₃⁴ + 2.0x₁x₂²x₃² + x₁x₂⁴ + 2.0x₁²x₃³ + 2.0x₁²x₂x₃² + 2.0x₁²x₂²x₃ + 2.0x₁²x₂³ + 2.0x₁³x₃² + 2.0x₁³x₂² + x₁⁴x₃ + x₁⁴x₂ + x₁⁵ + 1.0e-5x₃⁶ + 3.0e-5x₂²x₃⁴ + 3.0e-5x₂⁴x₃² + 1.0e-5x₂⁶ + 3.0e-5x₁²x₃⁴ + 6.0e-5x₁²x₂²x₃² + 3.0e-5x₁²x₂⁴ + 3.0e-5x₁⁴x₃² + 3.0e-5x₁⁴x₂² + 1.0e-5x₁⁶
Objective was scaled by the prefactor 1.0 + 2.0x₃² + 2.0x₂² + 2.0x₁² + x₃⁴ + 2.0x₂²x₃² + x₂⁴ + 2.0x₁²x₃² + 2.0x₁²x₂² + x₁⁴
1 equality constraint
1: -1.0 + x₁x₂x₃ = 0
3 nonnegative constraints
2: x₁ ≥ 0
3: x₂ ≥ 0
4: x₃ ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 3 atol = 3e-4
            end
        end
    end
end

@testset "POP 23 (USSR Olimpian 1989)" begin
    DynamicPolynomials.@polyvar x[1:3]
    Random.seed!(1234)
    prob = poly_problem((x[1] + x[2])*(x[2] + x[3]), nonneg=x, zero=[prod(x)*sum(x)-1], noncompact=(1e-5, 5))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: 1.0e-5 + 7.0e-5x₃² + x₂x₃ + 1.00007x₂² + x₁x₃ + x₁x₂ + 7.0e-5x₁² + 0.00021x₃⁴ + 5.0x₂x₃³ + 5.00042x₂²x₃² + 5.0x₂³x₃ + 5.00021x₂⁴ + 5.0x₁x₃³ + 5.0x₁x₂x₃² + 5.0x₁x₂²x₃ + 5.0x₁x₂³ + 0.00042x₁²x₃² + 5.0x₁²x₂x₃ + 5.00042x₁²x₂² + 5.0x₁³x₃ + 5.0x₁³x₂ + 0.00021x₁⁴ + 0.00035x₃⁶ + 10.0x₂x₃⁵ + 10.00105x₂²x₃⁴ + 20.0x₂³x₃³ + 20.00105x₂⁴x₃² + 10.0x₂⁵x₃ + 10.00035x₂⁶ + 10.0x₁x₃⁵ + 10.0x₁x₂x₃⁴ + 20.0x₁x₂²x₃³ + 20.0x₁x₂³x₃² + 10.0x₁x₂⁴x₃ + 10.0x₁x₂⁵ + 0.00105x₁²x₃⁴ + 20.0x₁²x₂x₃³ + 20.0021x₁²x₂²x₃² + 20.0x₁²x₂³x₃ + 20.00105x₁²x₂⁴ + 20.0x₁³x₃³ + 20.0x₁³x₂x₃² + 20.0x₁³x₂²x₃ + 20.0x₁³x₂³ + 0.00105x₁⁴x₃² + 10.0x₁⁴x₂x₃ + 10.00105x₁⁴x₂² + 10.0x₁⁵x₃ + 10.0x₁⁵x₂ + 0.00035x₁⁶ + 0.00035x₃⁸ + 10.0x₂x₃⁷ + 10.0014x₂²x₃⁶ + 30.0x₂³x₃⁵ + 30.0021x₂⁴x₃⁴ + 30.0x₂⁵x₃³ + 30.0014x₂⁶x₃² + 10.0x₂⁷x₃ + 10.00035x₂⁸ + 10.0x₁x₃⁷ + 10.0x₁x₂x₃⁶ + 30.0x₁x₂²x₃⁵ + 30.0x₁x₂³x₃⁴ + 30.0x₁x₂⁴x₃³ + 30.0x₁x₂⁵x₃² + 10.0x₁x₂⁶x₃ + 10.0x₁x₂⁷ + 0.0014x₁²x₃⁶ + 30.0x₁²x₂x₃⁵ + 30.0042x₁²x₂²x₃⁴ + 60.0x₁²x₂³x₃³ + 60.0042x₁²x₂⁴x₃² + 30.0x₁²x₂⁵x₃ + 30.0014x₁²x₂⁶ + 30.0x₁³x₃⁵ + 30.0x₁³x₂x₃⁴ + 60.0x₁³x₂²x₃³ + 60.0x₁³x₂³x₃² + 30.0x₁³x₂⁴x₃ + 30.0x₁³x₂⁵ + 0.0021x₁⁴x₃⁴ + 30.0x₁⁴x₂x₃³ + 30.0042x₁⁴x₂²x₃² + 30.0x₁⁴x₂³x₃ + 30.0021x₁⁴x₂⁴ + 30.0x₁⁵x₃³ + 30.0x₁⁵x₂x₃² + 30.0x₁⁵x₂²x₃ + 30.0x₁⁵x₂³ + 0.0014x₁⁶x₃² + 10.0x₁⁶x₂x₃ + 10.0014x₁⁶x₂² + 10.0x₁⁷x₃ + 10.0x₁⁷x₂ + 0.00035x₁⁸ + 0.00021x₃¹⁰ + 5.0x₂x₃⁹ + 5.00105x₂²x₃⁸ + 20.0x₂³x₃⁷ + 20.0021x₂⁴x₃⁶ + 30.0x₂⁵x₃⁵ + 30.0021x₂⁶x₃⁴ + 20.0x₂⁷x₃³ + 20.00105x₂⁸x₃² + 5.0x₂⁹x₃ + 5.00021x₂¹⁰ + 5.0x₁x₃⁹ + 5.0x₁x₂x₃⁸ + 20.0x₁x₂²x₃⁷ + 20.0x₁x₂³x₃⁶ + 30.0x₁x₂⁴x₃⁵ + 30.0x₁x₂⁵x₃⁴ + 20.0x₁x₂⁶x₃³ + 20.0x₁x₂⁷x₃² + 5.0x₁x₂⁸x₃ + 5.0x₁x₂⁹ + 0.00105x₁²x₃⁸ + 20.0x₁²x₂x₃⁷ + 20.0042x₁²x₂²x₃⁶ + 60.0x₁²x₂³x₃⁵ + 60.0063x₁²x₂⁴x₃⁴ + 60.0x₁²x₂⁵x₃³ + 60.0042x₁²x₂⁶x₃² + 20.0x₁²x₂⁷x₃ + 20.00105x₁²x₂⁸ + 20.0x₁³x₃⁷ + 20.0x₁³x₂x₃⁶ + 60.0x₁³x₂²x₃⁵ + 60.0x₁³x₂³x₃⁴ + 60.0x₁³x₂⁴x₃³ + 60.0x₁³x₂⁵x₃² + 20.0x₁³x₂⁶x₃ + 20.0x₁³x₂⁷ + 0.0021x₁⁴x₃⁶ + 30.0x₁⁴x₂x₃⁵ + 30.0063x₁⁴x₂²x₃⁴ + 60.0x₁⁴x₂³x₃³ + 60.0063x₁⁴x₂⁴x₃² + 30.0x₁⁴x₂⁵x₃ + 30.0021x₁⁴x₂⁶ + 30.0x₁⁵x₃⁵ + 30.0x₁⁵x₂x₃⁴ + 60.0x₁⁵x₂²x₃³ + 60.0x₁⁵x₂³x₃² + 30.0x₁⁵x₂⁴x₃ + 30.0x₁⁵x₂⁵ + 0.0021x₁⁶x₃⁴ + 20.0x₁⁶x₂x₃³ + 20.0042x₁⁶x₂²x₃² + 20.0x₁⁶x₂³x₃ + 20.0021x₁⁶x₂⁴ + 20.0x₁⁷x₃³ + 20.0x₁⁷x₂x₃² + 20.0x₁⁷x₂²x₃ + 20.0x₁⁷x₂³ + 0.00105x₁⁸x₃² + 5.0x₁⁸x₂x₃ + 5.00105x₁⁸x₂² + 5.0x₁⁹x₃ + 5.0x₁⁹x₂ + 0.00021x₁¹⁰ + 7.0e-5x₃¹² + x₂x₃¹¹ + 1.00042x₂²x₃¹⁰ + 5.0x₂³x₃⁹ + 5.00105x₂⁴x₃⁸ + 10.0x₂⁵x₃⁷ + 10.0014x₂⁶x₃⁶ + 10.0x₂⁷x₃⁵ + 10.00105x₂⁸x₃⁴ + 5.0x₂⁹x₃³ + 5.00042x₂¹⁰x₃² + x₂¹¹x₃ + 1.00007x₂¹² + x₁x₃¹¹ + x₁x₂x₃¹⁰ + 5.0x₁x₂²x₃⁹ + 5.0x₁x₂³x₃⁸ + 10.0x₁x₂⁴x₃⁷ + 10.0x₁x₂⁵x₃⁶ + 10.0x₁x₂⁶x₃⁵ + 10.0x₁x₂⁷x₃⁴ + 5.0x₁x₂⁸x₃³ + 5.0x₁x₂⁹x₃² + x₁x₂¹⁰x₃ + x₁x₂¹¹ + 0.00042x₁²x₃¹⁰ + 5.0x₁²x₂x₃⁹ + 5.0021x₁²x₂²x₃⁸ + 20.0x₁²x₂³x₃⁷ + 20.0042x₁²x₂⁴x₃⁶ + 30.0x₁²x₂⁵x₃⁵ + 30.0042x₁²x₂⁶x₃⁴ + 20.0x₁²x₂⁷x₃³ + 20.0021x₁²x₂⁸x₃² + 5.0x₁²x₂⁹x₃ + 5.00042x₁²x₂¹⁰ + 5.0x₁³x₃⁹ + 5.0x₁³x₂x₃⁸ + 20.0x₁³x₂²x₃⁷ + 20.0x₁³x₂³x₃⁶ + 30.0x₁³x₂⁴x₃⁵ + 30.0x₁³x₂⁵x₃⁴ + 20.0x₁³x₂⁶x₃³ + 20.0x₁³x₂⁷x₃² + 5.0x₁³x₂⁸x₃ + 5.0x₁³x₂⁹ + 0.00105x₁⁴x₃⁸ + 10.0x₁⁴x₂x₃⁷ + 10.0042x₁⁴x₂²x₃⁶ + 30.0x₁⁴x₂³x₃⁵ + 30.0063x₁⁴x₂⁴x₃⁴ + 30.0x₁⁴x₂⁵x₃³ + 30.0042x₁⁴x₂⁶x₃² + 10.0x₁⁴x₂⁷x₃ + 10.00105x₁⁴x₂⁸ + 10.0x₁⁵x₃⁷ + 10.0x₁⁵x₂x₃⁶ + 30.0x₁⁵x₂²x₃⁵ + 30.0x₁⁵x₂³x₃⁴ + 30.0x₁⁵x₂⁴x₃³ + 30.0x₁⁵x₂⁵x₃² + 10.0x₁⁵x₂⁶x₃ + 10.0x₁⁵x₂⁷ + 0.0014x₁⁶x₃⁶ + 10.0x₁⁶x₂x₃⁵ + 10.0042x₁⁶x₂²x₃⁴ + 20.0x₁⁶x₂³x₃³ + 20.0042x₁⁶x₂⁴x₃² + 10.0x₁⁶x₂⁵x₃ + 10.0014x₁⁶x₂⁶ + 10.0x₁⁷x₃⁵ + 10.0x₁⁷x₂x₃⁴ + 20.0x₁⁷x₂²x₃³ + 20.0x₁⁷x₂³x₃² + 10.0x₁⁷x₂⁴x₃ + 10.0x₁⁷x₂⁵ + 0.00105x₁⁸x₃⁴ + 5.0x₁⁸x₂x₃³ + 5.0021x₁⁸x₂²x₃² + 5.0x₁⁸x₂³x₃ + 5.00105x₁⁸x₂⁴ + 5.0x₁⁹x₃³ + 5.0x₁⁹x₂x₃² + 5.0x₁⁹x₂²x₃ + 5.0x₁⁹x₂³ + 0.00042x₁¹⁰x₃² + x₁¹⁰x₂x₃ + 1.00042x₁¹⁰x₂² + x₁¹¹x₃ + x₁¹¹x₂ + 7.0e-5x₁¹² + 1.0e-5x₃¹⁴ + 7.0e-5x₂²x₃¹² + 0.00021x₂⁴x₃¹⁰ + 0.00035x₂⁶x₃⁸ + 0.00035x₂⁸x₃⁶ + 0.00021x₂¹⁰x₃⁴ + 7.0e-5x₂¹²x₃² + 1.0e-5x₂¹⁴ + 7.0e-5x₁²x₃¹² + 0.00042x₁²x₂²x₃¹⁰ + 0.00105x₁²x₂⁴x₃⁸ + 0.0014x₁²x₂⁶x₃⁶ + 0.00105x₁²x₂⁸x₃⁴ + 0.00042x₁²x₂¹⁰x₃² + 7.0e-5x₁²x₂¹² + 0.00021x₁⁴x₃¹⁰ + 0.00105x₁⁴x₂²x₃⁸ + 0.0021x₁⁴x₂⁴x₃⁶ + 0.0021x₁⁴x₂⁶x₃⁴ + 0.00105x₁⁴x₂⁸x₃² + 0.00021x₁⁴x₂¹⁰ + 0.00035x₁⁶x₃⁸ + 0.0014x₁⁶x₂²x₃⁶ + 0.0021x₁⁶x₂⁴x₃⁴ + 0.0014x₁⁶x₂⁶x₃² + 0.00035x₁⁶x₂⁸ + 0.00035x₁⁸x₃⁶ + 0.00105x₁⁸x₂²x₃⁴ + 0.00105x₁⁸x₂⁴x₃² + 0.00035x₁⁸x₂⁶ + 0.00021x₁¹⁰x₃⁴ + 0.00042x₁¹⁰x₂²x₃² + 0.00021x₁¹⁰x₂⁴ + 7.0e-5x₁¹²x₃² + 7.0e-5x₁¹²x₂² + 1.0e-5x₁¹⁴
Objective was scaled by the prefactor 1.0 + 5.0x₃² + 5.0x₂² + 5.0x₁² + 10.0x₃⁴ + 20.0x₂²x₃² + 10.0x₂⁴ + 20.0x₁²x₃² + 20.0x₁²x₂² + 10.0x₁⁴ + 10.0x₃⁶ + 30.0x₂²x₃⁴ + 30.0x₂⁴x₃² + 10.0x₂⁶ + 30.0x₁²x₃⁴ + 60.0x₁²x₂²x₃² + 30.0x₁²x₂⁴ + 30.0x₁⁴x₃² + 30.0x₁⁴x₂² + 10.0x₁⁶ + 5.0x₃⁸ + 20.0x₂²x₃⁶ + 30.0x₂⁴x₃⁴ + 20.0x₂⁶x₃² + 5.0x₂⁸ + 20.0x₁²x₃⁶ + 60.0x₁²x₂²x₃⁴ + 60.0x₁²x₂⁴x₃² + 20.0x₁²x₂⁶ + 30.0x₁⁴x₃⁴ + 60.0x₁⁴x₂²x₃² + 30.0x₁⁴x₂⁴ + 20.0x₁⁶x₃² + 20.0x₁⁶x₂² + 5.0x₁⁸ + x₃¹⁰ + 5.0x₂²x₃⁸ + 10.0x₂⁴x₃⁶ + 10.0x₂⁶x₃⁴ + 5.0x₂⁸x₃² + x₂¹⁰ + 5.0x₁²x₃⁸ + 20.0x₁²x₂²x₃⁶ + 30.0x₁²x₂⁴x₃⁴ + 20.0x₁²x₂⁶x₃² + 5.0x₁²x₂⁸ + 10.0x₁⁴x₃⁶ + 30.0x₁⁴x₂²x₃⁴ + 30.0x₁⁴x₂⁴x₃² + 10.0x₁⁴x₂⁶ + 10.0x₁⁶x₃⁴ + 20.0x₁⁶x₂²x₃² + 10.0x₁⁶x₂⁴ + 5.0x₁⁸x₃² + 5.0x₁⁸x₂² + x₁¹⁰
1 equality constraint
1: -1.0 + x₁x₂x₃² + x₁x₂²x₃ + x₁²x₂x₃ = 0
3 nonnegative constraints
2: x₁ ≥ 0
3: x₂ ≥ 0
4: x₃ ≥ 0"
    if optimize
        :MosekSOS ∈ solvers && @test poly_optimize(:MosekSOS, prob).objective ≈ 2 atol = 1e-5
        :COPTMoment ∈ solvers && @test poly_optimize(:COPTMoment, prob).objective ≈ 2 atol = 1e-4
        #:HypatiaMoment ∈ solvers && @test poly_optimize(:HypatiaMoment, prob).objective ≈ 2 atol = 2e-4 # 80s
        # rest slow
    end
end

@testset "POP 24 (IMO 1990)" begin
    DynamicPolynomials.@polyvar x[1:4]
    Random.seed!(1234)
    prob = poly_problem(sum(x[i]*prod(sum(x[j] for j in 1:4 if j != k) for k in 1:4 if k != i) for i in 1:4) -
        prod(sum(x[i] for i in 1:4 if i != j) for j in 1:4)/3, nonneg=x, zero=[x[1]*x[2]+x[2]*x[3]+x[3]*x[4]+x[4]*x[1]-1],
        noncompact=(1e-5, 0))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 4 variables
Objective: 1.0e-5 + 3.0e-5x₄² + 3.0e-5x₃² + 3.0e-5x₂² + 3.0e-5x₁² + 1.00003x₄⁴ + 1.66667x₃x₄³ + 1.33339x₃²x₄² + 1.66667x₃³x₄ + 1.00003x₃⁴ + 1.66667x₂x₄³ + 3.66667x₂x₃x₄² + 3.66667x₂x₃²x₄ + 1.66667x₂x₃³ + 1.33339x₂²x₄² + 3.66667x₂²x₃x₄ + 1.33339x₂²x₃² + 1.66667x₂³x₄ + 1.66667x₂³x₃ + 1.00003x₂⁴ + 1.66667x₁x₄³ + 3.66667x₁x₃x₄² + 3.66667x₁x₃²x₄ + 1.66667x₁x₃³ + 3.66667x₁x₂x₄² + 5.0x₁x₂x₃x₄ + 3.66667x₁x₂x₃² + 3.66667x₁x₂²x₄ + 3.66667x₁x₂²x₃ + 1.66667x₁x₂³ + 1.33339x₁²x₄² + 3.66667x₁²x₃x₄ + 1.33339x₁²x₃² + 3.66667x₁²x₂x₄ + 3.66667x₁²x₂x₃ + 1.33339x₁²x₂² + 1.66667x₁³x₄ + 1.66667x₁³x₃ + 1.66667x₁³x₂ + 1.00003x₁⁴ + 1.0e-5x₄⁶ + 3.0e-5x₃²x₄⁴ + 3.0e-5x₃⁴x₄² + 1.0e-5x₃⁶ + 3.0e-5x₂²x₄⁴ + 6.0e-5x₂²x₃²x₄² + 3.0e-5x₂²x₃⁴ + 3.0e-5x₂⁴x₄² + 3.0e-5x₂⁴x₃² + 1.0e-5x₂⁶ + 3.0e-5x₁²x₄⁴ + 6.0e-5x₁²x₃²x₄² + 3.0e-5x₁²x₃⁴ + 6.0e-5x₁²x₂²x₄² + 6.0e-5x₁²x₂²x₃² + 3.0e-5x₁²x₂⁴ + 3.0e-5x₁⁴x₄² + 3.0e-5x₁⁴x₃² + 3.0e-5x₁⁴x₂² + 1.0e-5x₁⁶
1 equality constraint
1: -1.0 + x₃x₄ + x₂x₃ + x₁x₄ + x₁x₂ = 0
4 nonnegative constraints
2: x₁ ≥ 0
3: x₂ ≥ 0
4: x₃ ≥ 0
5: x₄ ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ 5.0625837 atol = solver==:SCSMoment ? 1e-3 : 2e-7
            end
        end
    end
end

@testset "POP 25 (IMO 2000)" begin
    DynamicPolynomials.@polyvar x[1:3]
    Random.seed!(1234)
    prob = poly_problem(-prod(x[mod(i-1, 1:3)]*x[i] - x[i] +1 for i in 1:3), nonneg=x, zero=[prod(x)-1], noncompact=(1e-5, 1))
    map_coefficients!(x -> round(x, digits=5), prob.objective)
    @test strRep(prob) == "Real-valued polynomial optimization problem in 3 variables
Objective: -0.99999 + x₃ + x₂ + x₁ - 0.99995x₃² - 2.0x₂x₃ - 0.99995x₂² - 2.0x₁x₃ - 2.0x₁x₂ - 0.99995x₁² + x₃³ + x₂x₃² + 2.0x₂²x₃ + x₂³ + 2.0x₁x₃² + 4.0x₁x₂x₃ + x₁x₂² + x₁²x₃ + 2.0x₁²x₂ + x₁³ + 0.0001x₃⁴ - 2.0x₂x₃³ + 0.0002x₂²x₃² - 2.0x₂³x₃ + 0.0001x₂⁴ - 2.0x₁x₃³ - 4.0x₁x₂x₃² - 4.0x₁x₂²x₃ - 2.0x₁x₂³ + 0.0002x₁²x₃² - 4.0x₁²x₂x₃ + 0.0002x₁²x₂² - 2.0x₁³x₃ - 2.0x₁³x₂ + 0.0001x₁⁴ + x₂²x₃³ + x₂⁴x₃ + x₁x₃⁴ + 4.0x₁x₂x₃³ + 2.0x₁x₂²x₃² + 4.0x₁x₂³x₃ + 2.0x₁²x₂x₃² + 2.0x₁²x₂²x₃ + x₁²x₂³ + x₁³x₃² + 4.0x₁³x₂x₃ + x₁⁴x₂ + 0.0001x₃⁶ + 0.0003x₂²x₃⁴ + 0.0003x₂⁴x₃² + 0.0001x₂⁶ - 2.0x₁x₂x₃⁴ - 2.0x₁x₂²x₃³ - 2.0x₁x₂³x₃² - 2.0x₁x₂⁴x₃ + 0.0003x₁²x₃⁴ - 2.0x₁²x₂x₃³ - 0.9994x₁²x₂²x₃² - 2.0x₁²x₂³x₃ + 0.0003x₁²x₂⁴ - 2.0x₁³x₂x₃² - 2.0x₁³x₂²x₃ + 0.0003x₁⁴x₃² - 2.0x₁⁴x₂x₃ + 0.0003x₁⁴x₂² + 0.0001x₁⁶ + x₁x₂²x₃⁴ + x₁x₂⁴x₃² + x₁²x₂x₃⁴ + x₁²x₂²x₃³ + x₁²x₂³x₃² + x₁²x₂⁴x₃ + x₁³x₂²x₃² + x₁⁴x₂x₃² + x₁⁴x₂²x₃ + 5.0e-5x₃⁸ + 0.0002x₂²x₃⁶ + 0.0003x₂⁴x₃⁴ + 0.0002x₂⁶x₃² + 5.0e-5x₂⁸ + 0.0002x₁²x₃⁶ - 0.9994x₁²x₂²x₃⁴ - 0.9994x₁²x₂⁴x₃² + 0.0002x₁²x₂⁶ + 0.0003x₁⁴x₃⁴ - 0.9994x₁⁴x₂²x₃² + 0.0003x₁⁴x₂⁴ + 0.0002x₁⁶x₃² + 0.0002x₁⁶x₂² + 5.0e-5x₁⁸ + 1.0e-5x₃¹⁰ + 5.0e-5x₂²x₃⁸ + 0.0001x₂⁴x₃⁶ + 0.0001x₂⁶x₃⁴ + 5.0e-5x₂⁸x₃² + 1.0e-5x₂¹⁰ + 5.0e-5x₁²x₃⁸ + 0.0002x₁²x₂²x₃⁶ + 0.0003x₁²x₂⁴x₃⁴ + 0.0002x₁²x₂⁶x₃² + 5.0e-5x₁²x₂⁸ + 0.0001x₁⁴x₃⁶ + 0.0003x₁⁴x₂²x₃⁴ + 0.0003x₁⁴x₂⁴x₃² + 0.0001x₁⁴x₂⁶ + 0.0001x₁⁶x₃⁴ + 0.0002x₁⁶x₂²x₃² + 0.0001x₁⁶x₂⁴ + 5.0e-5x₁⁸x₃² + 5.0e-5x₁⁸x₂² + 1.0e-5x₁¹⁰
Objective was scaled by the prefactor 1.0 + x₃² + x₂² + x₁²
1 equality constraint
1: -1.0 + x₁x₂x₃ = 0
3 nonnegative constraints
2: x₁ ≥ 0
3: x₂ ≥ 0
4: x₃ ≥ 0"
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ -0.997439 atol = solver==:SCSMoment ? 2e-4 : 1e-6
            end
        end
    end
end
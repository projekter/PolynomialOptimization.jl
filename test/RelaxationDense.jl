include("./shared.jl")

@testset "Some real-valued examples (scalar)" begin
    DynamicPolynomials.@polyvar x[1:2]
    prob = poly_problem(-(x[1] - 1)^2 - (x[1] - x[2])^2 - (x[2] - 3)^2,
                        nonneg=[1 - (x[1] - 1)^2, 1 - (x[1] - x[2])^2, 1 - (x[2] - 3)^2], perturbation=1e-3)
    if optimize
        for solver in solvers
            for (i, sol) in ((1, -3.), (2, -2.))
                @testset let i=i, solver=solver
                    @test poly_optimize(solver, prob, i).objective ≈ sol atol = 2e-2
                end
            end
        end
    end
end

@testset "Some real-valued examples (matrix)" begin
    DynamicPolynomials.@polyvar x[1:2]
    prob = poly_problem(-x[1]^2 - x[2]^2, psd=[[1-4x[1]*x[2] x[1]; x[1] 4-x[1]^2-x[2]^2]], perturbation=1e-4)
    if optimize
        for solver in solvers
            for i in 1:2
                @testset let i=i, solver=solver
                    @test poly_optimize(solver, prob, i).objective ≈ -4. atol = 1e-3
                end
            end
        end
    end
end

@testset "Some real-valued examples (matrix and equality)" begin
    DynamicPolynomials.@polyvar x[1:2]
    prob = poly_problem(-x[1]^2 - x[2]^2, zero=[x[1] + x[2] - 1], psd=[[1-4x[1]*x[2] x[1]; x[1] 4-x[1]^2-x[2]^2]])
    if optimize
        for solver in solvers
            for (i, sol) in ((1, -4.), (2, -3.904891578336841))
                @testset let i=i, solver=solver
                    res = poly_optimize(solver, prob, i)
                    @test res.objective ≈ sol atol = 5e-3
                    if i == 2 && solver != :SCSMoment
                        for extraction in (solver == :SpecBMSOS ? (:heuristic,) : (:mvhankel, :heuristic))
                            pas = poly_all_solutions(extraction, res)
                            @test isone(length(pas))
                            @test pas[1][1] ≈ [-0.8047780612688199, 1.804778061268820] atol = (solver === :SpecBMSOS ? 1e-4 : 1e-6)
                        end
                    end
                end
            end
        end
    end
end

@testset "Some complex-valued examples (equality)" begin
    DynamicPolynomials.@complex_polyvar z
    prob = poly_problem(z + conj(z), zero=[z * conj(z) - 1])
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob).objective ≈ -2. atol = (solver === :SpecBMSOS ? 1e-4 : 1e-7)
            end
        end
    end
end

@testset "Some complex-valued examples (inequalities)" begin
    DynamicPolynomials.@complex_polyvar z[1:2]
    prob = poly_problem(3 - z[1] * conj(z[1]) - 0.5im * z[1] * conj(z[2])^2 + 0.5im * z[2]^2 * conj(z[1]),
                        zero=[z[1] * conj(z[1]) - 0.25z[1]^2 - 0.25conj(z[1])^2 - 1,
                              z[1] * conj(z[1]) + z[2] * conj(z[2]) - 3,
                              1im * z[2] - 1im * conj(z[2])], nonneg=[z[2] + conj(z[2])])
    if optimize
        for solver in solvers
            @testset let solver=solver
                @test poly_optimize(solver, prob, 3).objective ≈ 0.42817465 atol = 1e-5 skip=solver===:SpecBMSOS
            end
        end
    end
end

@testset "Some complex-valued examples (matrix)" begin
    DynamicPolynomials.@complex_polyvar x[1:2]
    prob = poly_problem(-x[1] * conj(x[1]) - x[2] * conj(x[2]),
                        psd=[[1-2*(x[1]*x[2]+conj(x[1] * x[2])) x[1]; conj(x[1]) 4-x[1]*conj(x[1])-x[2]*conj(x[2])]],
                        perturbation=1e-4)
    if optimize
        for solver in solvers
            for i in 2:3
                @testset let i=i, solver=solver
                    @test poly_optimize(solver, prob, i).objective ≈ -4. atol = 5e-3 skip=solver===:SpecBMSOS
                end
            end
        end
    end
end
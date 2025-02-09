include("./shared.jl")
using LinearAlgebra

roundfn(tol) = x -> round(x, digits=tol)
hs(x) = x' * x
almost_equal(x::P, y::P; tol=6) where {P<:AbstractPolynomial} = iszero(map_coefficients!(roundfn(tol), x - y))
almost_equal(x::Vector{<:SimplePolynomial}, y::AbstractPolynomial; tol=6) =
    almost_equal(sum(let v=variables(y); p -> hs(PolynomialOptimization.change_backend(p, v)) end, x, init=zero(y)), y; tol)

@testset "No constraints" begin
    DynamicPolynomials.@polyvar x[1:3];
    DynamicPolynomials.@complex_polyvar z[1:2];

    for (str, obj, vars, success) in (
        ("Real-valued PSD", 1 + x[1]^4 + x[2]^4 + x[3]^4 + x[1]^2*x[2]^2 + x[1]^2*x[3]^2 + x[2]^2*x[3]^2 + x[2]*x[3], x, true),
        ("Complex-valued PSD", hs(1 + z[1]*z[2]) + hs(1 + z[1]^2) + hs(3 - z[2]), z, false),
        ("Real-valued quadratic", 9 + x[1] + x[1]^2, [x[1]], true),
        ("Complex-valued quadratic", 9 + z[1] + conj(z[1]) + z[1]*conj(z[1]), [z[1]], true),
        ("Scalar", one(polynomial_type(x)), typeof(x[1])[], true)
    )
        @testset "$str" begin
            prob = poly_problem(obj)
            rel = Relaxation.Dense(prob)
            g = PolynomialOptimization.change_backend(Relaxation.groupings(rel).obj[1], vars)
            for solver in solvers, representation in (RepresentationPSD(), RepresentationDD(), RepresentationSDD())
                @testset let solver=solver, representation=representation
                    res = solver === :SCSMoment ? poly_optimize(solver, rel; eps_abs=1e-8, eps_rel=1e-8, representation) :
                                                  poly_optimize(solver, rel; representation)
                    @test issuccess(res) == !(representation isa RepresentationDD) || success
                    if success
                        cert = SOSCertificate(res)
                        sosc = dot(g, cert.data[1][1], g)
                        @test almost_equal(cert[:objective, 1], sosc)
                        @test almost_equal(sosc, obj - res.objective)
                    end
                end
            end
        end
    end
end

# the extraction logic is the same for everything, so let us now just check the decomposition once

@testset "Nonnegative constraint" begin
    DynamicPolynomials.@polyvar x[1:2];
    obj = -(x[1] - 1)^2 - (x[1] - x[2])^2 - (x[2] - 3)^2
    nonneg = [1 - (x[1] - 1)^2, 1 - (x[1] - x[2])^2, 1 - (x[2] - 3)^2]
    prob = poly_problem(obj; nonneg)

    rel = Relaxation.Dense(prob, 2)
    groupings = Relaxation.groupings(rel)
    g₁ = PolynomialOptimization.change_backend(groupings.obj[1], x)
    gₙ = PolynomialOptimization.change_backend.(first.(groupings.nonnegs), (x,))

    for solver in solvers
        @testset let solver=solver
            res = poly_optimize(solver, rel)
            cert = SOSCertificate(res)
            sosc1 = dot(g₁, cert.data[1][1], g₁)
            @test almost_equal(cert[:objective, 1], sosc1)
            sosc2 = dot(gₙ[1], cert.data[2][1], gₙ[1])
            @test almost_equal(cert[:nonneg, 1, 1], sosc2)
            sosc3 = dot(gₙ[2], cert.data[3][1], gₙ[2])
            @test almost_equal(cert[:nonneg, 2, 1], sosc3)
            sosc4 = dot(gₙ[3], cert.data[4][1], gₙ[3])
            @test almost_equal(cert[:nonneg, 3, 1], sosc4)
            @test almost_equal(sosc1 + nonneg[1] * sosc2 + nonneg[2] * sosc3 + nonneg[3] * sosc4, obj - res.objective,
                tol=solver === :SCSMoment ? 2 : 6)
        end
    end
end

@testset "PSD constraint" begin
    DynamicPolynomials.@complex_polyvar z[1:2]
    obj = -z[1] * conj(z[1]) - z[2] * conj(z[2])
    psd = [1-2*(z[1]*z[2]+conj(z[1] * z[2])) z[1]; conj(z[1]) 4-z[1]*conj(z[1])-z[2]*conj(z[2])]
    prob = poly_problem(obj, psd=[psd])

    rel = Relaxation.Dense(prob, 3)
    groupings = Relaxation.groupings(rel)
    g₁ = PolynomialOptimization.change_backend(groupings.obj[1], z)
    g₂ = PolynomialOptimization.change_backend(groupings.psds[1][1], z)

    for solver in solvers
        @testset let solver=solver
            if solver == :COPTMoment
                res = poly_optimize(solver, rel, parameters=((COPT.COPT_DBLPARAM_FEASTOL, 1e-8),
                                                             (COPT.COPT_DBLPARAM_DUALTOL, 1e-8)))
            else
                res = poly_optimize(solver, rel)
            end
            cert = SOSCertificate(res)
            sosc1 = dot(g₁, cert.data[1][1], g₁)
            @test almost_equal(cert[:objective, 1], sosc1)
            ref = reshape(cert.data[2][1], (size(psd, 1), length(g₂), size(psd, 2), length(g₂)))
            # now do einsum("a, iajb, b", conj(g₂), cm, g₂). Note that the second system must come first (col major).
            sosc2 = zeros(polynomial_type(eltype(psd), Complex{Float64}), size(psd)...)
            for b in 1:length(g₂), j in axes(psd, 2), a in 1:length(g₂), i in axes(psd, 1)
                sosc2[i, j] += conj(g₂[a]) * ref[i, a, j, b] * g₂[b]
            end
            cmp = hs(PolynomialOptimization.change_backend.(cert[:psd, 1, 1], (z,)))
            @test size(cmp) == size(psd)
            @test all(splat(almost_equal), zip(cmp, sosc2))
            @test almost_equal(sosc1 + sum(conj(x) * y for (x, y) in zip(sosc2, psd)), obj - res.objective,
                tol=solver === :SCSMoment ? 3 : 6)
        end
    end
end

@testset "Equality constraint" begin
    DynamicPolynomials.@polyvar x[1:2]
    DynamicPolynomials.@complex_polyvar z[1:2]
    for (str, obj, zero, nonneg, deg, vars) in (
        ("Real-valued", -x[1]^2 - x[2]^2, [x[1] + x[2] - 1], [(1 - 4x[1]*x[2])*(4 - x[1]^2 - x[2]^2) - x[1]^2,
                                                              1 - 4x[1]*x[2], 4 - x[1]^2 - x[2]^2], 2, x),
        ("Complex-valued", 3 - z[1] * conj(z[1]) - 0.5im * z[1] * conj(z[2])^2 + 0.5im * z[2]^2 * conj(z[1]),
        [z[1] * conj(z[1]) - 0.25z[1]^2 - 0.25conj(z[1])^2 - 1, z[1] * conj(z[1]) + z[2] * conj(z[2]) - 3,
         1im * z[2] - 1im * conj(z[2])], [z[2] + conj(z[2])], 3, z)
    )
        @testset "$str" begin
            prob = poly_problem(obj; zero, nonneg)

            rel = Relaxation.Dense(prob, deg)
            groupings = Relaxation.groupings(rel)
            g₁ = PolynomialOptimization.change_backend(groupings.obj[1], vars)

            for solver in solvers
                @testset let solver=solver
                    res = solver === :SCSMoment ? poly_optimize(solver, rel, eps_abs=1e-8, eps_rel=1e-8) :
                                                  poly_optimize(solver, rel)
                    cert = SOSCertificate(res)
                    # Checking with a reference implementation for the equality constraints is not possible easily - the
                    # prefactors would include anything from the matrix g' * g; however, duplicates were filtered and the order
                    # is basically random (dependent on how the sets arranges it). This is why cert[:zero] actually does a lot
                    # of work, and we wouldn't have a chance but to do the same work, which is pointless for the test.
                    # But we can check the result.
                    if !isreal(prob) # not implemented
                        @test_throws ErrorException cert[:zero, 1, 1]
                    else
                        @test(almost_equal(sum(hs, PolynomialOptimization.change_backend.(cert[:objective, 1], (vars,))) +
                                           sum(zero[i] * PolynomialOptimization.change_backend(cert[:zero, i, 1], vars)
                                               for i in 1:length(zero)) +
                                           sum(nonneg[i] *
                                               sum(hs, PolynomialOptimization.change_backend.(cert[:nonneg, i, 1], (vars,)))
                                               for i in 1:length(nonneg)),
                                           obj - res.objective),
                              broken=solver===:COPTMoment) # TODO: fix COPT. But is it our fault or COPT's? The data just seems
                                                           # to be wrong...
                    end
                end
            end
        end
    end
end
function lancelot_solve(problem::PolynomialOptimization.RealProblem; feastol::Number=1e-5, gradtol::Number=1e-5,
    maxit::Integer=1000, verbose=false)
    @verbose_info("Compiling polynomials for fast evaluation")
    setup_static = @elapsed begin
        nvars = nvariables(problem)
        static_objective = StaticPolynomials.Polynomial(problem.objective)
        static_scalar_constrs = Vector{StaticPolynomials.Polynomial}(undef, length(problem.constr_zero) +
            length(problem.constr_nonneg))
        i = 1
        @inbounds for constr in problem.constr_zero
            static_scalar_constrs[i] = StaticPolynomials.Polynomial(constr)
            i += 1
        end
        @inbounds for constr in problem.constr_nonneg
            rmul!(coefficients(constr), -1)
            static_scalar_constrs[i] = StaticPolynomials.Polynomial(constr)
            rmul!(coefficients(constr), -1)
            i += 1
        end
        num_scalar = length(static_scalar_constrs)
        # In this implementation, we use the smallest eigenvalue (which is always computed numerically using LAPACK's spexv)
        # as the single nonnegative constraint. We can also calculate the gradient with moderate effort; however, it is only in
        # the case of a nondegenerate lowest subspace that this is actually _the_ gradient and not just a subgradient. So it
        # would be good to study the effects of this subtlety.
        # An alternative way would be to add n inequality constraints for an n-dimensional PSD matrix and calculate the
        # coefficients of the characteristic polynomial (LaBudde's method would probably be the best - the reduction to the
        # tridigonal form is done any also by spevx). However, we then have to calculate the derivatives of the coefficients,
        # and the question is how to do this efficiently (and not by reducing n submatrices).
        # We could also try to explicitly construct the coefficients in symbolic form as StaticPolynomials (which would be the
        # most efficient way to go for very small matrices), but for larger matrices, this will - unless the matrix is very
        # sparse - probably be the worst choice due to the high degrees.
        static_psd_constrs = Vector{StaticPolynomials.PolynomialSystem}(undef, length(problem.constr_psd))
        i = 1
        @inbounds for constr in problem.constr_psd
            sd = LinearAlgebra.checksquare(constr)
            static_psd_constrs[i] = StaticPolynomials.PolynomialSystem(
                trttp!('U', constr, Vector{eltype(constr)}(undef, packedsize(sd))), variables=variables(problem.objective)
            )
            i += 1
        end
    end
    @verbose_info("Compilation finished in ", setup_static, " seconds.")
    psd_buffer_dict = Dict{Int,SPMatrix{Float64,Vector{Float64},:U}}()
    psd_buffer = SPMatrix{Float64,Vector{Float64},:U}[
        let k::Int=size(constr, 1)
            get!(() -> SPMatrix{Float64}(undef, k), psd_buffer_dict, k)
        end for constr in problem.constr_psd
    ]
    jac_buffer_dict = Dict{Int,Matrix{Float64}}()
    jac_buffer = Matrix{Float64}[
        get!(() -> Matrix{Float64}(undef, StaticPolynomials.npolynomials(constr), nvars), jac_buffer_dict,
            StaticPolynomials.npolynomials(constr))
        for constr in static_psd_constrs
    ]
    function eval_fun end
    let static_objective=static_objective, num_scalar=num_scalar, static_scalar_constrs=static_scalar_constrs,
        static_psd_constrs=static_psd_constrs, psd_buffer=psd_buffer
        eval_fun(x) = StaticPolynomials.evaluate(static_objective, x)
        function eval_fun(x, i)
            @inbounds if i ≤ num_scalar
                return StaticPolynomials.evaluate(static_scalar_constrs[i], x)
            else
                i -= num_scalar
                StaticPolynomials.evaluate!(vec(psd_buffer[i]), static_psd_constrs[i], x)
                return -first(eigvals!(psd_buffer[i], 1:1))
            end
        end
    end
    function eval_grad end
    let static_objective=static_objective, num_scalar=num_scalar, static_scalar_constrs=static_scalar_constrs,
        static_psd_constrs=static_psd_constrs, psd_buffer=psd_buffer, jac_buffer=jac_buffer
        eval_grad(g, x) = StaticPolynomials.gradient!(g, static_objective, x)
        function eval_grad(g, x, i)
            @inbounds if i ≤ num_scalar
                StaticPolynomials.gradient!(g, static_scalar_constrs[i], x)
            else
                i -= num_scalar
                # ∇λₘᵢₙ(Gᵢ(x)) = λₘᵢₙ'(Gᵢ(x)) Gᵢ'(x)
                # now ∇Gᵢ is simple, and ∂λₘᵢₙ(S₀ + ϵ S)/∂ϵ = ⟨v₀, S v₀⟩, where vᵢ: min eigvec of S₀
                # derivative of smallest eigenvalue of S₀ + ϵ S: eigenvalues of ⟨v₀, S v₀⟩, where v₀: min eigvec of S₀
                # λₘᵢₙ'(Y) = [∂λₘᵢₙ(Y + ϵ Eⱼ)/∂ϵ = ⟨y₀, E y₀⟩]ⱼ where Eⱼ are the symmetric elementary matrices
                psd = psd_buffer[i]
                StaticPolynomials.evaluate_and_jacobian!(vec(psd), jac_buffer[i], static_psd_constrs[i], x)
                v₀ = @view(eigen!(psd, 1:1).vectors[:, 1]) # eigen! overwrites psd with eigenvalues, but allocates for vectors
                dλ = vec(psd) # so we can just reuse psd here
                k = 1
                for j in 1:size(psd, 2) # we calculate ⟨v₀, Eᵢⱼ v₀⟩ directly, using the particular form of Eᵢⱼ
                    for i in 1:j-1
                        dλ[k] = 2v₀[j] * v₀[i]
                        k += 1
                    end
                    dλ[k] = v₀[j]^2
                    k += 1
                end
                mul!(reshape(g, (1, length(g))), reshape(dλ, (1, length(dλ))), jac_buffer[i], -1, false)
            end
        end
    end
    # we only provide the Hessian if there are no PSD constraints (as we either have to give all or none)
    if isempty(static_psd_constrs)
        function eval_hess end
        let static_objective=static_objective, num_scalar=num_scalar, static_scalar_constrs=static_scalar_constrs
            eval_hess(H, x) = StaticPolynomials.hessian!(H, static_objective, x)
            function eval_hess(H, x, i)
                @assert(i ≤ num_scalar)
                return StaticPolynomials.hessian!(H, static_scalar_constrs[i], x)
            end
        end
    else
        eval_hess = missing
    end
    opti_local = let nvars=nvars, MY_FUN=eval_fun, MY_GRAD=eval_grad, MY_HESS=eval_hess, neq=length(problem.constr_zero),
        nin=length(problem.constr_nonneg) + length(problem.constr_psd), feastol=feastol, gradtol=gradtol,
        print_level=verbose ? 1 : 0, verbose=verbose, maxit=maxit
        (x0) -> begin
            f, it, status = LANCELOT_simple(nvars, x0, MY_FUN; MY_GRAD, MY_HESS, neq, nin, feastol, gradtol, maxit,
                print_level)
            verbose && @printf("NLP: %d, iterations: %d, best value: %.8g\n", status, it, f)
            return iszero(status) ? f : eltype(x0)(Inf), x0
        end
    end
    return opti_local
end

function Solver.poly_optimize(::Val{:LANCELOT}, problem::PolynomialOptimization.Problem; verbose::Bool=false, feastol=1e-5,
    gradtol=1e-5, maxit::Integer=1000)
    return lancelot_solve(problem; feastol, gradtol, maxit, verbose)
end
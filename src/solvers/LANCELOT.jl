function lancelot_solve(problem; feastol::Number=1e-5, gradtol::Number=1e-5, maxit::Integer=1000, verbose=false,
    print_level::Int=verbose ? 1 : 0)
    @assert(!problem.complex)
    @verbose_info("Compiling polynomials for fast evaluation")
    setup_static = @elapsed begin
        nvars = length(problem.variables)
        static_objective = StaticPolynomials.Polynomial(problem.objective, variables=problem.variables)
        # Important: this relies on the fact that in the definition of PolyConstraintType, the equalities come before the
        # inequalities, which is as LANCELOT expects it. Additionally, we use that PSD constraints come last
        gettype = x -> x.type
        constraints = sort(problem.constraints, by=gettype)
        if isempty(constraints)
            equality_constraints = 0
            inequality_constraints = 0
        else
            # as we specify by to access the type field, we must also make our search element such that they have a type field.
            cmp = first(constraints)
            equality_constraints = searchsortedfirst(constraints, typeof(cmp)(pctNonneg, cmp.constraint, cmp.basis),
                by=gettype) -1
            inequality_constraints = searchsortedfirst(@view(constraints[equality_constraints+1:end]),
                typeof(cmp)(pctPSD, cmp.constraint, cmp.basis), by=gettype) -1
        end
        scalar_constraints = equality_constraints + inequality_constraints
        psd_constraints = length(constraints) - scalar_constraints
        static_scalar_constrs = StaticPolynomials.Polynomial[
                                    StaticPolynomials.Polynomial(c.type == pctNonneg ? -c.constraint : c.constraint,
                                        variables=problem.variables) # we could just make everything negative, but let's keep
                                                                     # close to the original problem
                                    for c in @view(constraints[begin:scalar_constraints])
                                ]
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
        static_psd_constr = StaticPolynomials.PolynomialSystem[let sd = LinearAlgebra.checksquare(c.constraint)
                                StaticPolynomials.PolynomialSystem(
                                    trttp!('U', c.constraint, Vector{eltype(c.constraint)}(undef, (sd * (sd +1)) >> 1))
                                )
                            end for c in @view(constraints[scalar_constraints+1:end])]
    end
    @verbose_info("Compilation finished in ", setup_static, " seconds.")
    psd_buffer_dict = Dict{Int,PackedMatrix{Float64,Vector{Float64}}}()
    psd_buffer = PackedMatrix{Float64,Vector{Float64},:U}[let k::Int=size(c.constraint, 1)
                                            get!(() -> PackedMatrix{Float64}(undef, k), psd_buffer_dict, k)
                                       end for c in @view(constraints[scalar_constraints+1:end])]
    jac_buffer_dict = Dict{Int,Matrix{Float64}}()
    jac_buffer = Matrix{Float64}[get!(() -> Matrix{Float64}(undef, StaticPolynomials.npolynomials(c), nvars), jac_buffer_dict,
                                      StaticPolynomials.npolynomials(c)) for c in static_psd_constr]
    function eval_fun end
    let static_objective=static_objective, scalar_constraints=scalar_constraints, static_scalar_constrs=static_scalar_constrs,
        static_psd_constr=static_psd_constr, psd_buffer=psd_buffer
        eval_fun(x) = StaticPolynomials.evaluate(static_objective, x)
        function eval_fun(x, i)
            @inbounds if i ≤ scalar_constraints
                return StaticPolynomials.evaluate(static_scalar_constrs[i], x)
            else
                i -= scalar_constraints
                StaticPolynomials.evaluate!(vec(psd_buffer[i]), static_psd_constr[i], x)
                return -first(eigvals!(psd_buffer[i], 1:1))
            end
        end
    end
    function eval_grad end
    let static_objective=static_objective, scalar_constraints=scalar_constraints, static_scalar_constrs=static_scalar_constrs,
        static_psd_constr=static_psd_constr, psd_buffer=psd_buffer, jac_buffer=jac_buffer
        eval_grad(g, x) = StaticPolynomials.gradient!(g, static_objective, x)
        function eval_grad(g, x, i)
            @inbounds if i ≤ scalar_constraints
                StaticPolynomials.gradient!(g, static_scalar_constrs[i], x)
            else
                i -= scalar_constraints
                # ∇λₘᵢₙ(Gᵢ(x)) = λₘᵢₙ'(Gᵢ(x)) Gᵢ'(x)
                # now ∇Gᵢ is simple, and ∂λₘᵢₙ(S₀ + ϵ S)/∂ϵ = ⟨v₀, S v₀⟩, where vᵢ: min eigvec of S₀
                # derivative of smallest eigenvalue of S₀ + ϵ S: eigenvalues of ⟨v₀, S v₀⟩, where v₀: min eigvec of S₀
                # λₘᵢₙ'(Y) = [∂λₘᵢₙ(Y + ϵ Eⱼ)/∂ϵ = ⟨y₀, E y₀⟩]ⱼ where Eⱼ are the symmetric elementary matrices
                psd = psd_buffer[i]
                StaticPolynomials.evaluate_and_jacobian!(vec(psd), jac_buffer[i], static_psd_constr[i], x)
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
    if iszero(psd_constraints)
        function eval_hess end
        let static_objective=static_objective, scalar_constraints=scalar_constraints,
            static_scalar_constrs=static_scalar_constrs
            eval_hess(H, x) = StaticPolynomials.hessian!(H, static_objective, x)
            function eval_hess(H, x, i)
                @assert(i ≤ scalar_constraints)
                return StaticPolynomials.hessian!(H, static_scalar_constrs[i], x)
            end
        end
    else
        eval_hess = missing
    end
    opti_local = let nvars=nvars, MY_FUN=eval_fun, MY_GRAD=eval_grad, MY_HESS=eval_hess, neq=equality_constraints,
        nin=inequality_constraints + psd_constraints, feastol=feastol, gradtol=gradtol, print_level=print_level,
        verbose=verbose, maxit=maxit
        (x0) -> begin
            f, it, status = LANCELOT_simple(nvars, x0, MY_FUN; MY_GRAD, MY_HESS, neq, nin, feastol, gradtol, maxit,
                print_level)
            verbose && @printf("NLP: %d, iterations: %d, best value: %.8g\n", status, it, f)
            return iszero(status) ? f : eltype(x0)(Inf), x0
        end
    end
    return opti_local
end

"""
    sparse_optimize(:LANCELOT, state; verbose, feastol, gradtol)

Alias for [`poly_optimize`](@ref), regardless of the sparsity pattern.
"""
function sparse_optimize(::Val{:LANCELOT}, state::AbstractSparsity; verbose::Bool=false, feastol=1e-5, gradtol=1e-5,
    maxit::Integer=1000)
    problem = sparse_problem(state)
    return lancelot_solve(problem; feastol, gradtol, maxit, verbose)
end
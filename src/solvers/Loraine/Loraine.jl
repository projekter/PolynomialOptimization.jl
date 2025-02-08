# See the original code under MIT license at https://github.com/kocvara/Loraine.jl; however, the solver is so badly implemented
# that we bascially build it again from scratch. This is still not fully optimal as there are some temporary allocations in the
# iterations that should be avoided, but it is so much more efficient than the original code. (Plus, it actually works on huge
# problems, as it doesn't try to allocate matrices that are never needed.)
# Eventually, this should go into a pull request; then, we'll delete this module and put it just as an interface to the weakref
# as all the other solvers. However, at the moment, all features that are not needed were removed: direct solver,
# rank-1 data matrices, detailed timing. Data is always sparse.
module Loraine

export solve!

using LinearAlgebra, SparseArrays, Printf, PositiveFactorizations, ConjugateGradients
using ..Solvers: EfficientCholmod
using SuiteSparse: CHOLMOD

include("./Enums.jl")
include("./DataTypes.jl")
include("./Preconditioners.jl")
include("./PredictorCorrector.jl")

"""
    solve!(solver)

Solves a freshly initialized solver. Note that calling this function twice will produce the same results (unless parameters are
changed), as the initial state is always restored at the beginning of the call.
"""
function solve!(solver::Solver)
    if solver.verb > VERBOSITY_NONE
        println(" *** IP STARTS")
        if solver.verb < VERBOSITY_FULL
            println(" it        obj         error     cg_iter")
        else
            println(" it        obj         error      err1      err2      err3      err4      err5      err6    cg_pre  cg_cor")
        end
    end

    initial_point!(solver)

    orig_aa = solver.aamat
    sqrtn = sqrt(solver.model.n)
    if solver.preconditioner == PRECONDITIONER_NONE
        preconditioner = PreconditionerNone()
    elseif solver.preconditioner == PRECONDITIONER_HALPHA
        preconditioner = PreconditionerAlpha(solver)
    else
        preconditioner = PreconditionerBeta(solver)
    end

    tottime = @elapsed while solver.status == STATUS_UNKNOWN
        ip_step!(solver, preconditioner)

        solver.tol_cg = max(solver.tol_cg * solver.tol_cg_up, solver.tol_cg_min)

        check_convergence!(solver)

        if solver.preconditioner == PRECONDITIONER_HYBRID && preconditioner isa PreconditionerBeta
            if (10solver.cg_iter_cor > solver.erank * solver.model.nlmi * sqrtn &&
                60solver.iter > sqrtn) || solver.cg_iter_cor > 100
                preconditioner = PreconditionerAlpha(solver)
                solver.aamat = AMAT_IDENTITY
                if solver.verb > VERBOSITY_NONE
                    println("Switching to preconditioner H_alpha")
                end
            end
        end
    end

    solver.aamat = orig_aa
    if solver.verb > VERBOSITY_NONE
        @printf(" *** Total CG iterations: %8.0d \n", solver.cg_iter_tot)
        solver.status == STATUS_OPTIMAL && @printf(" *** Optimal solution found in %8.2f seconds\n", tottime)
    end
end

function initial_point!(solver::Solver{T}) where {T<:Real}
    if solver.initpoint != INITPOINT_LORAINE
        nrmb2 = sqrt(sum(x -> (one(T) + abs(x))^2, solver.model.b, init=zero(T)))
    end
    fill!(solver.y, zero(T))
    @inbounds for i in 1:solver.model.nlmi
        if solver.initpoint == INITPOINT_LORAINE
            copyto!(solver.X[i], I)
            η = solver.model.n
        else
            f = nrmb2 / (one(T) + norm(solver.model.A[i]))
            sqrtm = sqrt(solver.model.coneDims[i])
            copyto!(solver.X[i], sqrt(solver.model.coneDims[i]) * max(one(T), sqrtm * f) * I)
            η = max(sqrtm, one(T) + max(f, norm(solver.model.C[i])))
        end
        copyto!(solver.S[i], η * I)
    end

    dd = solver.model.nlin
    if dd > 0
        if solver.initpoint == INITPOINT_LORAINE
            fill!(solver.X_lin, one(T))
            fill!(solver.S_lin, one(T))
        else
            fill!(solver.X_lin, max(one(T),
                                    maximum((one(T) + abs(bᵢ)) / (one(T) + norm(A_linᵢ))
                                            for (bᵢ, A_linᵢ) in zip(solver.model.b, eachrow(solver.model.A_lin)), init=typemin(T))::T))
            fill!(solver.S_lin, max(one(T),
                                    max(maximum(norm, eachrow(solver.model.A_lin), init=typemin(T))::T,
                                        norm(solver.model.c_lin)) / sqrt(dd)))
        end
        fill!(solver.S_lin_inv, inv(first(solver.S_lin)))
    end

    solver.sigma = T(3)
    solver.DIMACS_error = one(T)
    solver.cg_iter_tot = 0
    solver.iter = 0
    solver.status = STATUS_UNKNOWN

    fill!.(solver.RNT, zero(T))
    fill!(solver.RNT_lin, zero(T))

    return solver
end

function ip_step!(solver::Solver{T,I}, @nospecialize(preconditioner::PreconditionerUnion{T,I})) where {T<:Real,I<:Integer}
    solver.iter += 1
    if solver.iter > solver.maxit
        solver.status = STATUS_ITERATION_LIMIT
        solver.verb > VERBOSITY_NONE &&
            println("WARNING: Stopped by iteration limit (stopping status = STATUS_ITERATION_LIMIT)")
        return solver
    end
    solver.cg_iter_pre = 0
    solver.cg_iter_cor = 0

    find_mu!(solver)
    prepare_W!(solver)
    predictor!(solver, preconditioner)
    sigma_update!(solver)
    corrector!(solver, preconditioner)

    return solver
end

function find_mu!(solver::Solver{T}) where {T<:Real}
    mu = sum(t -> dot(Symmetric(t[1], :U), Symmetric(t[2], :U)), zip(solver.X, solver.S), init=zero(T)) +
         dot(solver.X_lin, solver.S_lin)
    solver.mu = mu / (sum(solver.model.coneDims, init=0) + solver.model.nlin)
    return solver.mu
end

function prepare_W!(solver::Solver{T}) where {T<:Real}
    @inbounds for i in 1:solver.model.nlmi
        Ctmp = cholesky(Positive, solver.X[i])
        copyto!(solver.Si[i], solver.S[i])
        CtmpS = cholesky!(Positive, solver.Si[i])
        _, Dtmp, V = svd!(mul!(solver.W[i], CtmpS.L', Ctmp.L)) # we overwrite W[i] later on anyway

        copyto!(solver.D[i], Dtmp)
        Di2 = solver.DDsi[i] # just as a temporary
        try
            Di2 .= inv.(sqrt.(Dtmp))
        catch
            println("WARNING: Numerical difficulties, giving up")
            solver.status = STATUS_NUMERICAL
            mul!(solver.G[i], Ctmp.L, V)
        else
            rmul!(mul!(solver.G[i], Ctmp.L, V), Diagonal(Di2))
        end

        ldiv!(factorize(solver.G[i]), copyto!(solver.Gi[i], I))
        mul!(solver.W[i], solver.G[i], solver.G[i]')
        LinearAlgebra.inv!(CtmpS) # puts it in solver.Si[i]
        DDtmp = Ctmp.factors
        mul!(DDtmp, solver.S[i], solver.G[i])
        try
            solver.DDsi[i] .= inv.(sqrt.(dot.(eachcol(solver.G[i]), eachcol(DDtmp))))
        catch err
            println("WARNING: Numerical difficulties, giving up")
            fill!(solver.DDsi[i], one(T))
            solver.status = STATUS_NUMERICAL
            return
        end
    end
    solver.Si_lin .= inv.(solver.S_lin)

    return
end

function sigma_update!(solver::Solver{T}) where {T<:Real}
    step_pred = min(minimum(solver.alpha, init=typemax(T)), solver.alpha_lin,
                    minimum(solver.beta, init=typemax(T)), solver.beta_lin)
    if solver.mu > T(1//1_000_000)
        expon_used = step_pred .< inv(sqrt(T(3))) ? one(T) : max(solver.expon, T(3) * step_pred^2)
    else
        expon_used = max(one(T), min(solver.expon, T(3) * step_pred^2))
    end
    tmp1 = sum(t -> dot(Symmetric(t[1], :U), Symmetric(t[2], :U)), zip(solver.Xn, solver.Sn), init=zero(T))
    if tmp1 < 0
        solver.sigma = T(8//10)
    else
        tmp2 = dot(solver.Xn_lin, solver.Sn_lin)
        tmp12 = (tmp1 + tmp2) / (sum(solver.model.coneDims, init=0) + solver.model.nlin)
        solver.sigma = min(one(T), (tmp12 / solver.mu) ^ expon_used)
    end

    return solver.sigma
end

function check_convergence!(solver::Solver{T}) where {T<:Real}
    # DIMACS error evaluation
    solver.err1 = norm(solver.Rp) / (one(T) + norm(solver.model.b))
    solver.err2, solver.err3, solver.err4, solver.err6 = zero(T), zero(T), zero(T), zero(T)
    CX = sum(t -> dot(Symmetric(t[1], :U), Symmetric(t[2], :U)), zip(solver.model.C, solver.X), init=zero(T))
    by = dot(solver.model.b, solver.y)
    @inbounds for i in 1:solver.model.nlmi
        Xᵢ, Sᵢ, Cᵢ = Symmetric(solver.X[i], :U), Symmetric(solver.S[i], :U), Symmetric(solver.model.C[i], :U)
        solver.err2 += max(zero(T), -eigmin(Xᵢ,) / (one(T) + norm(solver.model.b)))
        solver.err3 += norm(Symmetric(solver.Rd[i], :U)) / (one(T) + norm(Cᵢ))
        solver.err4 += max(zero(T), -eigmin(Sᵢ) / (one(T) + norm(Cᵢ)))
        # err5 += (dot(Cᵢ, Xᵢ) - by) / (one(T) + abs(dot(Cᵢ, Xᵢ)) + abs(by)
        solver.err6 += dot(Sᵢ, Xᵢ) / (one(T) + abs(dot(Cᵢ, Xᵢ)) + abs(by))
    end
    if solver.model.nlin > 0
        dX = dot(solver.model.c_lin, solver.X_lin)
        solver.err2 += max(zero(T), -minimum(solver.X_lin, init=typemax(T)) / (one(T) + norm(solver.model.b)))
        solver.err3 += norm(solver.Rd_lin) / (one(T) + norm(solver.model.c_lin))
        solver.err4 += max(zero(T), -minimum(solver.S_lin, init=typemax(T)) / (one(T) + norm(solver.model.c_lin)))
        solver.err5 = (CX + dX - by) / (one(T) + abs(CX) + abs(by))
        solver.err6 += dot(solver.S_lin , solver.X_lin) / (one(T) + abs(dX) + abs(by))
    else
        solver.err5 = (CX - by) / (one(T) + abs(CX) + abs(by))
    end

    solver.DIMACS_error = solver.err2 + solver.err3 + solver.err4 + abs(solver.err5) + solver.err6
    if solver.model.nlmi > 0
        solver.DIMACS_error += solver.err1
    end
    if solver.verb > VERBOSITY_NONE && solver.status == STATUS_UNKNOWN
        if solver.verb > VERBOSITY_SHORT
            @printf("%3.0d %16.8e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %7.0d %7.0d\n", solver.iter,
                -by, solver.DIMACS_error, solver.err1, solver.err2, solver.err3, solver.err4, solver.err5, solver.err6,
                solver.cg_iter_pre, solver.cg_iter_cor)
        else
            @printf("%3.0d %16.8e %9.2e %9.0d\n", solver.iter, -by, solver.DIMACS_error,
                solver.cg_iter_pre + solver.cg_iter_cor)
        end
    end

    if solver.DIMACS_error < solver.eDIMACS
        solver.status = STATUS_OPTIMAL
        if solver.verb > VERBOSITY_NONE
            println("Primal objective: ", -by)
            println("Dual objective:   ", solver.model.nlin > 0 ? -CX - dX : -CX)
        end
    end

    if solver.DIMACS_error > T(1e25)
        solver.status = STATUS_INFEASIBLE
        solver.verb > VERBOSITY_NONE &&
            println("WARNING: Problem probably infeasible (stopping status = STATUS_INFEASIBLE)")
    elseif solver.DIMACS_error > T(1e25) || abs(by) > T(1e25)
        solver.status = STATUS_INFEASIBLE_OR_UNBOUNDED
        solver.verb > VERBOSITY_NONE &&
            println("WARNING: Problem probably unbounded or infeasible (stopping status = STATUS_INFEASIBLE_OR_UNBOUNDED)")
    end

    return
end

end
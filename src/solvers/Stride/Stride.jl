# This is an implementation of the STRIDE solver, https://doi.org/10.1007/s10107-022-01912-6
module Stride

using LinearAlgebra, SparseArrays, StandardPacked, Printf
using ...PolynomialOptimization: @assert, @inbounds, @verbose_info, poly_problem, poly_optimize, poly_solutions, Result,
    Relaxation
using ..Solvers: EfficientCholmod
using ...Solver: MomentVector
import Optim

export stride_solve

struct Data{R<:Real,F<:Factorization,Fmt,SS<:Union{SparseMatrixCSC{R,Int},UniformScaling{Bool}}}
    d::Int
    A::Vector{SparseMatrixCSC{R,Int}}
    Astar::Vector{SparseMatrixCSC{R,Int}}
    AAstarinv::F
    C::Vector{SPMatrix{R,SparseVector{R,Int},Fmt}}
    equality_subspace::Tuple{SS,SparseVector{R,Int}}
    b::SparseVector{R,Int}
    offset::R
    X::Vector{SPMatrix{R,Vector{R},Fmt}} # each dimensions nᵢ × nᵢ
    S::Vector{SPMatrix{R,Vector{R},Fmt}}
    tmpm1::Vector{SPMatrix{R,Vector{R},Fmt}}
    tmpm2::Vector{SPMatrix{R,Vector{R},Fmt}}
    y::Vector{R} # length d
    tmpv::Vector{R}
    tmpvlong::Vector{R} # length packedsize(max(nᵢ))
end

function psdproject!(m::SPMatrix{R}) where {R}
    munscaled = packed_unscale!(m)
    eigs = eigen(munscaled, zero(R), typemax(R))
    fill!(munscaled, zero(eltype(m)))
    for (eval, evec) in zip(eigs.values, eachcol(eigs.vectors))
        eval > zero(R) && spr!(eval, evec, munscaled)
    end
    return packed_scale!(munscaled)
end

mutable struct PolyLBFGSState{R,S<:Optim.LBFGSState} <: Optim.AbstractOptimizerState
    const state::S
    const ssd::Data{R}
    const tol::R
    success::Bool
end

Optim.reset!(method, state::PolyLBFGSState, obj, x) = Optim.reset!(method, state.state, obj, x)
Optim.update_state!(d, state::PolyLBFGSState, method::Optim.LBFGS) = Optim.update_state!(d, state.state, method)
Optim.update_h!(d, state::PolyLBFGSState, method::Optim.LBFGS) = Optim.update_h!(d, state.state, method)
Optim.trace!(tr, d, state::PolyLBFGSState, iteration, method::Optim.LBFGS, options, curr_time=time()) =
    Optim.trace!(tr, d, state.state, iteration, method, options, curr_time)
function Base.getproperty(state::PolyLBFGSState, prop::Symbol)
    if prop === :state || prop === :ssd || prop === :tol || prop === :success
        return getfield(state, prop)
    else
        return getproperty(getfield(state, :state), prop)
    end
end
function Base.setproperty!(state::PolyLBFGSState, prop::Symbol, v)
    if prop === :success
        return setfield!(state, prop, v)
    else
        return setfield!(getfield(state, :state), prop, v)
    end
end

function Optim.assess_convergence(state::PolyLBFGSState{R}, d, ::Optim.Options) where {R}
    ssd = state.ssd
    Astar, Z, W, ξ = ssd.Astar, ssd.X, ssd.S, state.state.x
    for (Wᵢ, Zᵢ, Astarᵢ) in zip(W, Z, Astar)
        copyto!(Wᵢ, Zᵢ)
        mul!(Wᵢ, Astarᵢ, ξ, -one(R), -one(R))
        psdproject!(Wᵢ)
    end
    # x_converged, f_converged, g_converged, f_increased - we disregard all this and use the criterion from Algorithm 6
    state.success = converged = ηproj(ssd, W, ξ) ≤ state.tol
    return converged, converged, converged, false
end

function stride_solve(ssd::Data{R}, relaxation::Relaxation.AbstractRelaxation, data; verbose::Bool=false,
    tol::R=1e-8, maxiter::Integer=10,
    tol_init::R=1e-4, maxiter_init::Integer=1000, verbose_init::Bool=false,
    tol_sgsapg::R=1e-12, maxiter_sgsapg::Integer=1000, verbose_sgsapg::Bool=false,
    tol_lbfgs::R=1e-12, mem_lbfgs::Integer=10, maxiter_lbfgs::Integer=1000, verbose_lbfgs::Bool=false,
    ls_lbfgs=Optim.LineSearches.MoreThuente(),
    σ::R=10.,
    verbose_local::Bool=false, kwargs_local::Dict=Dict(),
    opti_local=begin
        @verbose_info("No local optimization specified; automatically constructing one (disable by passing opti_local=nothing).")
        setup_static = @elapsed(ol = poly_optimize(Val(:LANCELOT), poly_problem(relaxation); verbose=verbose_local,
            kwargs_local...))
        @verbose_info("Local optimization constructed in ", setup_static, " seconds.")
        ol
    end,
    opti_round=nothing) where {R<:Real}
    (tol_init > 0 && maxiter_init > 0 && tol > 0 && σ > 0) || throw(ArgumentError("Invalid arguments"))
    # A: output, [X₁, ...], α, β -> α output + β A(X₁, ...)
    # Astar: output, y, i, α, β -> α output + β Xᵢ(y)
    # AAstarinv: output, y -> ̃y
    Astar, b, C, X, S, y, tmpm = ssd.Astar, ssd.b, ssd.C, ssd.X, ssd.S, ssd.y, ssd.tmpm1
    # Input: Generate
    # - an initial point (X₀, y₀, S₀) via spadmm (i.e., Algorithm 4). Not necessarily feasible.
    spadmm!(ssd, σ, (1 + sqrt(5)) / 2, tol_init, maxiter_init, verbose_init)
    all(M -> isposdef(M, tol), S) || error("An initial point was generated that is not positive semidefinite")
    # - a nondecreasing positive sequence {σₖ}; choose σₖ ≥ O(√k) to get convergence according to Theorem 2.
    #   > choose σₖ = sqrt(k) ... according to the IEEE paper, choose 10 (constant)
    # - a nonnegative sequence {ϵₖ} such that {kϵₖ} is summable,
    #   > looks like it is not required
    # - a stopping criterion η: 𝕊₊ⁿ × ℝᵐ × 𝕊₊ⁿ → ℝ₊ with a tolerance Tol > 0
    # - V = {X₀} if X₀ is feasible, V = ∅ otherwise. We don't store the Xs, but just their objectives
    #   > make it feasible by projecting onto the positive part as a first step.
    @verbose_info("Post-validating the initial point")
    for (Xᵢ, Sᵢ) in zip(X, S)
        psdproject!(Xᵢ)
        psdproject!(Sᵢ)
    end
    Vopt = dot(b, ssd.y)
    #   > now we have positive X and S, but feasibility also requires the linear constraints to be satisfied. Check it for the
    #     dual problem:
    ηd = zero(R)
    for (Astarᵢ, Sᵢ, Cᵢ, tmpmᵢ) in zip(Astar, S, C, tmpm)
        copyto!(tmpmᵢ, Sᵢ)
        mul!(tmpmᵢ, Astarᵢ, y, true, true)
        axpy!(-one(R), Cᵢ, tmpmᵢ)
        ηd += norm(tmpmᵢ)
    end
    if ηd > tol
        # no, we are infeasible. So Vopt is not an appropriate value for the last good result.
        @verbose_info("Initial candidate generated by ADMM+ is infeasible with error ", ηd)
        Vopt = R(Inf)
    else
        @verbose_info("Initial candidate generated by ADMM+ is feasible")
    end
    CX = Vopt
    # - a positive integer r ∈ [1, n] - this is the rank
    # - a positive constant ϵ > 0. According to the IEEE paper, choose 1e-12

    if !isnothing(opti_local)
        groupings = Relaxation.groupings(relaxation)
        maxfn = gr -> maximum(length, gr, init=0)
        tmpmon = Vector{R}(undef, max(maxfn(groupings.obj), maximum(maxfn, groupings.nonnegs, init=0),
            maximum(maxfn, groupings.psds, init=0)))
        sqrt2 = sqrt(R(2))
    end
    local lbfgs_fg!
    let Z=X, A=ssd.A, Astar=ssd.Astar, b=ssd.b
        # ϕ(ξ) := ‖Π(A* ξ + Z)‖^2/2 - ⟨b, ξ⟩
        # ∇ϕ(ξ) = AΠ(A* ξ + Z) - b
        lbfgs_fg! = Optim.only_fg!((F, G, ξ) -> begin
            ϕ = zero(R)
            isnothing(G) || (G .= .-b)
            for (i, (tmpmᵢ, Zᵢ, Aᵢ, Astarᵢ)) in enumerate(zip(tmpm, Z, A, Astar))
                copyto!(tmpmᵢ, Zᵢ)
                mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
                psdproject!(tmpmᵢ)
                isnothing(F) || (ϕ += norm(tmpmᵢ)^2)
                isnothing(G) || mul!(G, Aᵢ, tmpmᵢ, true, true)
            end
            isnothing(F) || return R(1/2)*ϕ - dot(b, ξ)
            nothing
        end)
    end
    lbfgs = Optim.LBFGS(m=mem_lbfgs, linesearch=ls_lbfgs)
    lbfgs_options = Optim.Options(g_tol=tol_lbfgs, show_trace=verbose_lbfgs, show_every=50, iterations=maxiter_lbfgs)
    lbfgs_d = Optim.OnceDifferentiable(lbfgs_fg!, y, zero(R))
    lbfgs_state = PolyLBFGSState(Optim.LBFGSState(
        y, # x
        similar(y), # x_previous
        similar(y), # g_previous
        Vector{R}(undef, mem_lbfgs), # rho
        [similar(y) for _ in 1:mem_lbfgs], # dx_history
        [similar(y) for _ in 1:mem_lbfgs], # dg_history
        fill(R(NaN), length(y)), # dx
        fill(R(NaN), length(y)), # # dg
        fill(R(NaN), length(y)), # u
        R(NaN), # f_x_previous
        similar(y), # twoloop_q
        Vector{R}(undef, mem_lbfgs), # twoloop_alpha
        0, # pseudo_iteration
        similar(y), # s
        similar(y), # x_ls
        one(R)
    ), ssd, tol_lbfgs, false)

    status = :max_iter
    σₖ = zero(R)
    # Iterate the following steps for k = 1, ...:
    @inbounds for k in 1:maxiter
        σₖ += σ
        # Step 1 (Projection). Compute (̄Xₖ, yₖ, Sₖ) satisfying (5) and (6), where ̄Xₖ is an inexact projection of Xₖ₋₁ - σₖ C onto
        # the primal feasible set of the SDP, by running algorithm sgsapg and algorithm mlbfgs, warmstarted by (yₖ₋₁, Sₖ₋₁).
        # → Section 4
        @inbounds for (Xᵢ, Cᵢ) in zip(X, C)
            axpy!(-σₖ, Cᵢ, Xᵢ)
        end
        sgsapg!(ssd, tol_sgsapg, maxiter_sgsapg, verbose_sgsapg)
        verbose_lbfgs && @verbose_info("Entering limited-memory BFGS method")
        # we have to duplicate the behavior of initial_state here
        Optim.reset!(lbfgs, lbfgs_state, lbfgs_d, y)
        copyto!(lbfgs_state.state.x_previous, y) # not strictly necessary (else it would be part of reset!)
        lbfgs_state.success = false
        Optim.optimize(lbfgs_d, y, lbfgs, lbfgs_options, lbfgs_state)
        if !lbfgs_state.success
            for (Wᵢ, Zᵢ, Astarᵢ) in zip(S, X, Astar)
                copyto!(Wᵢ, Zᵢ)
                mul!(Wᵢ, Astarᵢ, y, -one(R), -one(R))
                psdproject!(Wᵢ)
            end
            ηproj(ssd, S, y) # to calculate X (which is stored in tmpm2)
        end
        # else checking the termination criteria already constructed the state
        # We got an updated y and S, now we have to construct X = Π(A*ξ + Z), where Z = Xₖ₋₁ - σₖ C ≡ X

        for (Xᵢ, Astarᵢ) in zip(X, Astar)
            mul!(Xᵢ, Astarᵢ, y, true, true)
            psdproject!(Xᵢ)
        end

        # (̄Xₖ, yₖ, Sₖ) = (X(W, ξ), 1/σₖ ξ, 1/σₖ W) [note that ξ and y are aliased, as well as W and S]
        σinvₖ = inv(σₖ)
        rmul!(y, σinvₖ)
        rmul!.(S, σinvₖ)

        # Step 2 (Certification). If η(̄Xₖ, yₖ, Sₖ) < Tol, output (̄Xₖ, yₖ, Sₖ) and stop.
        CX, by, ηp, ηd, ηg = η(ssd, Val(true))
        (verbose_sgsapg || !isnothing(opti_local) || isone(k)) &&
            @verbose_info("Iteration | Primal objective | Dual objective | Primal infeasibility | Dual infeasibility | Relative duality gap")
        @verbose_info(@sprintf("%9d | %16g | %14g | %20g | %18g | %20g", k, CX, by, ηp, ηd, ηg))
        if max(ηp, ηg, ηd) ≤ tol
            @verbose_info("Successfully converged.")
            Vopt = CX
            status = :ok
            break
        end

        # Step 3 (Acceleration): Compute ̂Xₖ via rounding, local NLP optimization and lifting, starting from ̄Xₖ:
        if !isnothing(opti_local)
            # Original description:
            # ̄Xₖ = ∑ᵢ λᵢ vᵢ vᵢᵀ (truncate after r)
            # ̄xₖᵢ = rounding(vᵢ). Problem-dependent, but can first normalize vᵢ such that the constant component is 1, then
            #                     extract the order-1 entries. If easily possible, project to the feasible set.
            # This would be problematic, to the claim in the Yang paper that "the algorithms...can be extended...in a
            # straightforward way to solve sparse relaxations. We now have multiple X, contrary to the Yang paper. If those
            # just came from constraints, we would simply disregard all but the moment matrix. However, we allow for generic
            # sparsity structure, i.e., there may be multiple sub-moment matrices with possibly overlapping variables. To
            # account for this, we will use our solution extraction heuristic. This is different from the eigenvector
            # decomposition method!
            # (Also note that the poor-man's option to add a dense first-order moment matrix for each clique suggested
            # previously by Wang et al. [CS-TSSOS] appears to be rather pointless.)
            @verbose_info("Entering local optimization")
            res = Result(relaxation, :StrideMoment, NaN64, nothing, :interrupted, CX)
            res.moments = MomentVector(relaxation, data, X, nothing)
            solutions = isnothing(opti_round) ? poly_solutions(:heuristic, res; verbose) : opti_round(res; verbose)
            # ̂xₖᵢ = nlp(̄xₖᵢ)
            # ̂xₖ = argmin p(̂xₖᵢ) over all { ̂xₖᵢ : i = 1, ..., r }
            best = typemax(R)
            local best_solution
            for solution in solutions
                for i in eachindex(solution)
                    if isnan(solution[i])
                        solution[i] = rand()
                    end
                end
                val, pos = opti_local(solution) # note: we REQUIRE solutions returned by opti_local to be feasible (should
                                                # return (Inf, arbitrary) if not)!
                # opti_local is allowed to modify its parameter, e.g., by writing pos into solution.
                if val < best
                    best_solution = pos
                    best = val
                end
            end
            # ̂Xₖ = proj(monomial lifting of ̂xₖ)
            # Step 4 (Policy). Choose the better candidate in {̄Xₖ, ̂Xₖ} to update Xₖ:
            # Xₖ = ̂Xₖ if ⟨C, ̂Xₖ⟩ < min( ⟨C, ̄Xₖ⟩, min_{X ∈ V} { ⟨C, X⟩ : X ∈ V } ) - ϵ ∧ ̂Xₖ feasible
            # We decide to be a bit more careful here, which is necessary due to the absence of BFGS. Our SGSAPG result may be
            # very infeasible, so it is not fair to compare to this value unless both feasibilities are sufficienty well under
            # control.
            if ηp < 100tol && ηd < 100tol
                bound = min(Vopt, CX)
            else
                bound = Vopt
            end
            if best < bound - R(1//100_000_000_000)
                @verbose_info("Found better optimization candidate with objective ", best)
                i = 1
                @inbounds for grtype in ((groupings.obj,), groupings.nonnegs, groupings.psds),
                    groups::typeof(groupings.obj) in grtype, grᵢ in groups
                    tmpmonᵢ = @view(tmpmon[1:length(grᵢ)])
                    tmpmonᵢ .= (best_solution,) .|> grᵢ
                    Xᵢ = X[i]
                    idx = 1
                    for col in 1:length(grᵢ)
                        Xᵢ[idx] = tmpmonᵢ[col]^2
                        idx += 1
                        for row in col+1:length(grᵢ)
                            Xᵢ[idx] = sqrt2 * tmpmonᵢ[col] * tmpmonᵢ[row]
                            idx += 1
                        end
                    end
                    i += 1
                end
                # If Xₖ = ̂Xₖ, set V ← V ∪ { ̂Xₖ }
                Vopt = best
                CX = best
            else
                @verbose_info("Optimization did not improve")
            end
        end
        #    = ̄Xₖ otherwise (nop, we already did this)
    end
    @verbose_info("Optimization completed")

    return status, CX
end

function spadmm!(ssd::Data{R}, σ::R, γ::R, tol::R, maxiter::Integer, verbose::Bool) where {R<:Real}
    @assert(0 ≤ γ ≤ 2 && σ > 0 && tol ≥ 0 && maxiter > 0)
    A, Astar, AAstarinv, C, b, X, S, y = ssd.A, ssd.Astar, ssd.AAstarinv, ssd.C, ssd.b, ssd.X, ssd.S, ssd.y
    tmpm = ssd.tmpm1
    @verbose_info("Entering ADMM+ with penalty σ = ", σ, ", step length γ = ", γ, ", and tolerance ", tol,
        "\nIteration | Primal objective | Dual objective | Primal infeasibility | Dual infeasibility | Relative duality gap")
    # Initialization: Initial points X₀ = S₀ = 0 ∈ 𝕊ⁿ
    fill!.(X, zero(R))
    fill!.(S, zero(R))
    σ⁻¹ = inv(σ)
    γσ = γ * σ
    local y
    # Iterate the following steps for k = 1, ...
    @inbounds for k in 1:maxiter
        # Step 1: Compute
        #    ̂yₖ₊₁ = (A A*)⁻¹(b/σ - A(Xₖ/σ + Sₖ - C))
        for (i, (Xᵢ, Sᵢ, Cᵢ, Aᵢ)) in enumerate(zip(X, S, C, A))
            # We don't need Sᵢ any more, so let's overwrite it
            Sᵢ .+= Xᵢ .* σ⁻¹
            axpy!(-one(R), Cᵢ, Sᵢ) # we use axpy! for all sparse changes instead of broadcasting, as this can exploit the
                                   # sparsity pattern. In contrast, axpy! seems to be much worse than broadcasting for dense
                                   # vectors.
            mul!(y, Aᵢ, Sᵢ, -one(R), !isone(i))
        end
        axpy!(σ⁻¹, b, y)
        ldiv!(AAstarinv, y)
        # Step 2: Compute
        #    Sₖ₊₁ = ( Π(Xₖ + σ(A* ̂yₖ₊₁ - C)) - (Xₖ + σ(A* ̂yₖ₊₁ - C)) ) / σ
        #             ^ projection onto PSD, use LAPACK.syevx
        for (Xᵢ, Sᵢ, tmpmᵢ, Astarᵢ, Cᵢ) in zip(X, S, tmpm, Astar, C)
            mul!(tmpmᵢ, Astarᵢ, y, σ, false)
            axpy!(-σ, Cᵢ, tmpmᵢ)
            tmpmᵢ .+= Xᵢ
            copyto!(Sᵢ, tmpmᵢ)
            # We effectively project onto the (absolute of the) negative part, which is done by a positive projection, since
            # the positive part is expected to have low rank
            psdproject!(Sᵢ)
            Sᵢ .= (Sᵢ .- tmpmᵢ) .* σ⁻¹
        end
        # Step 3: Compute
        #    yₖ₊₁ = (A A*)⁻¹(b/σ - A(Xₖ/σ + Sₖ₊₁ - C))
        for (i, (Xᵢ, Sᵢ, Cᵢ, tmpmᵢ, Aᵢ)) in enumerate(zip(X, S, C, tmpm, A))
            # This time, we need the Sᵢ again (as well as the Xs), so we need the temporary storage
            tmpmᵢ .= Xᵢ .* σ⁻¹ .+ Sᵢ
            axpy!(-one(R), Cᵢ, tmpmᵢ)
            mul!(y, Aᵢ, tmpmᵢ, -one(R), !isone(i))
        end
        axpy!(σ⁻¹, b, y)
        ldiv!(AAstarinv, y)
        # Step 4: Compute
        #    Xₖ₊₁ = Xₖ + γ σ (Sₖ₊₁ + A* yₖ₊₁ - C)
        for (Xᵢ, Sᵢ, Cᵢ, Astarᵢ) in zip(X, S, C, Astar)
            Xᵢ .+= γσ .* Sᵢ
            axpy!(-γσ, Cᵢ, Xᵢ)
            mul!(Xᵢ, Astarᵢ, y, γσ, true)
        end
        # Until η(Xₖ₊₁, yₖ₊₁, Sₖ₊₁) ≤ tol (eq. 28)
        CX, by, ηp, ηd, ηg = η(ssd, Val(true))
        success = max(ηp, ηg, ηd) ≤ tol
        if verbose && (success || k == maxiter || isone(k % 20))
            @verbose_info(@sprintf("%9d | %16g | %14g | %20g | %18g | %20g", k, CX, by, ηp, ηd, ηg))
        end
        # Output Xₖ₊₁, yₖ₊₁, Sₖ₊₁
        if success
            @verbose_info("Found an initial point up to the given tolerance")
            return
        end
    end
    @verbose_info("Terminated initial point finding due to maxiter")
end

function sgsapg!(ssd::Data{R}, tol::R, maxiter::Integer, verbose::Bool) where {R<:Real}
    @assert(tol ≥ 0) # also all(W .⪰ 0), but this is too expensive to check
    A, Astar, AAstarinv, b, Z, W, ξ = ssd.A, ssd.Astar, ssd.AAstarinv, ssd.b, ssd.X, ssd.S, ssd.y
    tmpm1, tmpm2, tmpv = ssd.tmpm1, ssd.tmpm2, ssd.tmpv
    @verbose_info("Entering sGS-based accelerated proximal gradient method\nIteration |     residue")
    # Initialization: Set ̃W₁ = W₀ and t₁ = 1
    copyto!.(tmpm1, W)
    t = 1.
    # Iterate the following steps for k = 1, ...
    k = 1
    @inbounds while true
        # Step 1 (sGS update): Compute
        #    ̃ξₖ = (A A*)⁻¹(b - A(Z) - A(̃Wₖ))
        for (i, (Aᵢ, Zᵢ, tmpm1ᵢ)) in enumerate(zip(A, Z, tmpm1))
            tmpm1ᵢ .+= Zᵢ
            mul!(tmpv, Aᵢ, tmpm1ᵢ, -one(R), !isone(i))
        end
        tmpv .+= b
        copyto!.(tmpm1, W) # we need to back up the previous W
        ldiv!(ξ, AAstarinv, tmpv)
        #    Wₖ = Π(-A*̃ξₖ - Z) = Π(A*ξₖ + Z) - (A*ξₖ + Z)
        for (Zᵢ, Wᵢ, tmpm2ᵢ, Astarᵢ) in zip(Z, W, tmpm2, Astar)
            copyto!(Wᵢ, Zᵢ)
            mul!(Wᵢ, Astarᵢ, ξ, true, true)
            copyto!(tmpm2ᵢ, Wᵢ)
            psdproject!(Wᵢ)
            Wᵢ .-= tmpm2ᵢ
        end
        #    ξₖ = (A A*)⁻¹(b - A(Z) - A(Wₖ))
        for (i, (Aᵢ, Zᵢ, Wᵢ)) in enumerate(zip(A, Z, W))
            mul!(tmpv, Aᵢ, Zᵢ, -one(R), !isone(i))
            mul!(tmpv, Aᵢ, Wᵢ, -one(R), true)
        end
        tmpv .+= b
        ldiv!(ξ, AAstarinv, tmpv)

        # Until ηproj(Wₖ, ξₖ) ≤ tol
        ηp = ηproj(ssd, W, ξ) # ηproj will use tmpm2, tmpvlong, and tmpv as temporaries
        success = ηp ≤ tol
        if verbose && (success || k == maxiter || isone(k % 20))
            @verbose_info(@sprintf("%9d | %11g", k, ηp))
        end
        if success
            return # Output Wₖ, ξₖ
        elseif k == maxiter
            @verbose_info("Maximum iteration count reached")
            return
        end

        # Step 2: Compute
        #    tₖ₊₁ = (1 + sqrt(1 + 4tₖ^2))/2
        newt = (one(R) + sqrt(one(R) + R(4)*t^2)) / R(2)
        #    ̃Wₖ₊₁ = Wₖ + (tₖ-1)/tₖ₊₁ (Wₖ - Wₖ₋₁)
        frac = (t - one(R)) / newt
        t = newt
        # now Wₖ is W, Wₖ₋₁ is tmpm1, and we need to fill ̃Wₖ₊₁, which will be tmpm1
        for (Wᵢ, tmpm1ᵢ) in zip(W, tmpm1)
            tmpm1ᵢ .= (one(R) + frac) .* Wᵢ .- frac .* tmpm1ᵢ
        end

        k += 1
    end
end

function η(ssd::Data{R}, ::Val{details}=Val(false)) where {R,details}
    A, Astar, C, b, X, S, y, tmpm, tmpv = ssd.A, ssd.Astar, ssd.C, ssd.b, ssd.X, ssd.S, ssd.y, ssd.tmpm1, ssd.tmpv
    # here, note that b = (1, 0, ..., 0)
    # ηp(X) = ‖A(X) - b‖/(1 + ‖b‖)
    for (i, (Aᵢ, Xᵢ)) in enumerate(zip(A, X))
        mul!(tmpv, Aᵢ, Xᵢ, true, !isone(i))
    end
    @inbounds tmpv .-= b
    ηp = norm(tmpv) / (1 + norm(b))
    # ηd(y, S) = ‖A*(y) + S - C‖/(1 + ‖C‖)
    # ηg(X, y) = |⟨C, X⟩ - ⟨b, y⟩|/(1 + |⟨C, X⟩| + |⟨b, y⟩|)
    ηdnumsq, normCsq, CX, by = zero(R), zero(R), zero(R), dot(b, y)
    @inbounds for (Sᵢ, Xᵢ, Cᵢ, tmpmᵢ, Astarᵢ) in zip(S, X, C, tmpm, Astar)
        mul!(tmpmᵢ, Astarᵢ, y, true, false)
        tmpmᵢ .+= Sᵢ
        axpy!(-one(R), Cᵢ, tmpmᵢ)
        CX += dot(Cᵢ, Xᵢ)
        ηdnumsq += LinearAlgebra.norm(tmpmᵢ)^2
        normCsq += LinearAlgebra.norm(Cᵢ)^2
    end
    ηd = sqrt(ηdnumsq) / (1 + sqrt(normCsq))
    ηg = abs(CX - by) / (1 + abs(CX) + abs(by))
    return details ? (ssd.offset - CX, ssd.offset - by, ηp, ηd, ηg) : max(ηp, ηd, ηg)
end

function ηproj(ssd::Data{R,<:Factorization,Fmt}, W, ξ, recompute_AstarξplusZ::Bool=true) where {R,Fmt}
    A, Astar, b, Z, X = ssd.A, ssd.Astar, ssd.b, ssd.X, ssd.tmpm2
    tmpmvec, tmpv = ssd.tmpvlong, ssd.tmpv
    # X(W, ξ) := A*ξ + W + Z
    # ηproj := max( ‖A(X(W, ξ)) - b‖, ‖X(W, ξ) - Π(A*ξ + Z)‖ )
    norm2 = zero(R)
    if recompute_AstarξplusZ
        @inbounds for (Xᵢ, Wᵢ, Zᵢ, Astarᵢ) in zip(X, W, Z, Astar)
            tmpmᵢ = let s = size(Xᵢ, 1)
                SPMatrix(s, view(tmpmvec, 1:packedsize(s)), Fmt)
            end
            copyto!(tmpmᵢ, Zᵢ)
            mul!(tmpmᵢ, Astarᵢ, ξ, true, true)
            Xᵢ .= tmpmᵢ .+ Wᵢ
            psdproject!(tmpmᵢ)
            j = 0
            for q in 1:size(tmpmᵢ, 2)
                @simd for k in 1:q-1
                    norm2 += 2(Xᵢ[j+k] - tmpmᵢ[j+k])^2
                end
                norm2 += (Xᵢ[j+q] - tmpmᵢ[j+q])^2
                j += q
            end
        end
    else
        # if this is called from the BFGS function, the algorithm guarantees that norm2 is indeed exactly zero. But we still
        # need X.
        @inbounds for (Xᵢ, Wᵢ, Zᵢ, Astarᵢ) in zip(X, W, Z, Astar)
            copyto!(Xᵢ, Zᵢ)
            mul!(Xᵢ, Astarᵢ, ξ, true, true)
            Xᵢ .+= Wᵢ
        end
    end
    for (i, (Aᵢ, Xᵢ)) in enumerate(zip(A, X))
        mul!(tmpv, Aᵢ, Xᵢ, true, !isone(i))
    end
    axpy!(-one(R), b, tmpv)
    norm1 = norm(tmpv)
    return max(norm1, sqrt(norm2))
end

end
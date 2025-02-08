struct AOperator{T<:Real,I<:Integer}
    solver::Solver{T,I}
end

function (t::AOperator{T})(Ax::Vector{T}, x::Vector{T}) where {T<:Real}
    # Ax = ∑ᵢ Aᵢ * (Wᵢ * mat(x * Aᵢ) * Wᵢ) + A_lin * ((X_lin .* S_lin_inv) .* (x * A_lin))
    init = true
    @inbounds for i in 1:t.solver.model.nlmi
        Aᵢ, Wᵢ, dimᵢ, tmp₁, tmp₂ = t.solver.model.A[i], t.solver.W[i], t.solver.model.coneDims[i], t.solver.delX[i],
            t.solver.delS[i]
        tmp₁v = Base.ReshapedArray(tmp₁, (length(tmp₁),), ()) # faster than vec, which allocates a new Vector (although it
                                                              # aliases, but the struct itself is allocated)
        mul!(tmp₁v, Aᵢ', x)
        mul!(tmp₂, Wᵢ, tmp₁)
        mul!(tmp₁, tmp₂, Wᵢ)
        mul!(Ax, Aᵢ, tmp₁v, true, !init)
        init = false
    end
    if t.solver.model.nlin > 0
        tmp = t.solver.delX_lin
        mul!(tmp, t.solver.model.A_lin', x)
        tmp .*= t.solver.X_lin .* t.solver.S_lin_inv
        mul!(Ax, t.solver.model.A_lin, tmp, true, true)
    end
    return Ax
end

function predictor!(solver::Solver{T,I}, @nospecialize(preconditioner::PreconditionerUnion{T,I})) where {T<:Real,I<:Integer}
    # rₚ = b - Aᵀ X; R_d = C - S - ∑ᵢ yᵢ Aᵢ
    copyto!(solver.Rp, solver.model.b)
    @inbounds for i in 1:solver.model.nlmi
        mul!(solver.Rp, solver.model.A[i], vec(solver.X[i]), -one(T), one(T))
        Rdᵢ, Cᵢ, Sᵢ = solver.Rd[i], solver.model.C[i], solver.S[i]
        for j in 1:size(Rdᵢ, 2)
            @simd for k in 1:j
                Rdᵢ[j, k] = Rdᵢ[k, j] = Cᵢ[k, j] - Sᵢ[k, j]
            end
        end
        mul!(vec(Rdᵢ), solver.model.A[i]', solver.y, -one(T), one(T))
    end
    if solver.model.nlin > 0
        mul!(solver.Rp, solver.model.A_lin, solver.X_lin, -one(T), one(T))
        solver.Rd_lin .= solver.model.c_lin .- solver.S_lin
        mul!(solver.Rd_lin, solver.model.A_lin', solver.y, -one(T), one(T))
    end

    h = copyto!(solver.rhs, solver.Rp)  # RHS for the Hessian equation: r = rₚ + Aᵀ vec[ W (R_d + S) W ]
    @inbounds for i in 1:solver.model.nlmi
        tmp₁, tmp₂ = solver.delX[i], solver.delS[i]
        tmp₁ .= solver.Rd[i] .+ solver.S[i]
        mul!(tmp₂, tmp₁, solver.W[i])
        mul!(tmp₁, solver.W[i], tmp₂)
        mul!(h, solver.model.A[i], Base.ReshapedArray(tmp₁, (length(tmp₁),), ()), one(T), one(T))
    end
    if solver.model.nlin > 0
        mul!(h, solver.model.A_lin, solver.X_lin .* solver.Si_lin .* solver.Rd_lin .+ solver.X_lin, one(T), one(T))
    end

    # solving the linear system
    A = AOperator(solver)
    fill!(solver.dely, zero(T))
    prepare_prec!(solver, preconditioner)
    # ConjugateGradients.jl needs `tol` to be `Float64`,
    # maybe we can fix this in that package but in the mean time, we just convert the tolerance to `Float64`
    exit_code, num_iters = cg!(A, h, solver.dely; tol=Float64(solver.tol_cg), maxIter=10000, precon=preconditioner)

    solver.cg_iter_pre += num_iters
    solver.cg_iter_tot += num_iters

    find_step!(solver, true)
    return solver
end

function corrector!(solver::Solver{T,I}, @nospecialize(preconditioner::PreconditionerUnion{T,I})) where {T<:Real,I<:Integer}
    h = copyto!(solver.rhs, solver.Rp) # RHS for the linear system
    σμ = solver.sigma * solver.mu
    @inbounds for i in 1:solver.model.nlmi
        tmp = mul!(solver.delX[i], solver.Rd[i], solver.G[i])
        tmp2 = mul!(solver.delS[i], solver.G[i]', tmp)
        Dᵢ = solver.D[i]
        @simd for j in 1:length(Dᵢ)
            tmp2[j, j] += Dᵢ[j] - σμ / Dᵢ[j]
        end
        tmp2 .-= solver.RNT[i]
        mul!(tmp, tmp2, solver.G[i]')
        mul!(tmp2, solver.G[i], tmp)
        mul!(h, solver.model.A[i], vec(tmp2), true, true)
    end
    for j in 1:solver.model.nlin
        axpy!(solver.Si_lin[j] * (solver.X_lin[j] * solver.Rd_lin[j] + solver.delX_lin[j] * solver.delS_lin[j] - σμ) +
                solver.X_lin[j], @view(solver.model.A_lin[:, j]), h)
    end

    # solving the linear system
    A = AOperator(solver)
    fill!(solver.dely, zero(T))
    exit_code, num_iters = cg!(A, h, solver.dely; tol=Float64(solver.tol_cg), maxIter=10000, precon=preconditioner)

    solver.cg_iter_cor += num_iters
    solver.cg_iter_tot += num_iters

    find_step!(solver, false)
    return solver
end

function find_step!(solver::Solver{T}, predict::Bool) where {T}
    @inbounds for i in 1:solver.model.nlmi
        delSᵢ, delXᵢ, Wᵢ, Xᵢ, Gᵢ, DDsiᵢ = solver.delS[i], solver.delX[i], solver.W[i], solver.X[i], solver.G[i], solver.DDsi[i]
        copyto!(delSᵢ, solver.Rd[i])
        mul!(vec(delSᵢ), solver.model.A[i]', solver.dely, -one(T), true)
        # temporaries:
        # - predict: Xn and RNT are overwritten later on anyway, so we can use them as temporaries
        # - correct: Si and RNT are only required once in this block, then they can be used.
        if predict
            tmp = mul!(solver.Xn[i], delSᵢ, Wᵢ')
            copyto!(delXᵢ, Xᵢ)
            mul!(delXᵢ, Wᵢ, tmp, -one(T), -one(T))
        else
            delXᵢ .= (solver.sigma * solver.mu) .* solver.Si[i] .- Xᵢ
            tmp = mul!(solver.Si[i], delSᵢ, Wᵢ')
            mul!(delXᵢ, Wᵢ, tmp, -one(T), true)
            mul!(tmp, solver.RNT[i], solver.G[i]')
            mul!(delXᵢ, Gᵢ, tmp, true, true)
        end

        # determining steplength to stay feasible
        mul!(tmp, delXᵢ, solver.Gi[i]')
        XXX = mul!(solver.RNT[i], solver.Gi[i], tmp)
        XXX .= DDsiᵢ' .* XXX .* DDsiᵢ
        mimiX = eigvals!(Symmetric(XXX, :U), 1:1)[1]
        solver.alpha[i] = mimiX .> T(-1//1_000_000) ? T(99//100) : min(one(T), -solver.tau / mimiX)

        mul!(tmp, delSᵢ, Gᵢ)
        mul!(XXX, Gᵢ', tmp)
        XXX .= DDsiᵢ' .* XXX .* DDsiᵢ
        mimiS = eigvals!(Symmetric(XXX, :U), 1:1)[1]
        solver.beta[i] = mimiS > T(-1//1_000_000) ? T(99//100) : min(one(T), -solver.tau / mimiS)
    end

    if solver.model.nlin > 0
        find_step_lin!(solver, predict)
    else
        solver.alpha_lin = one(T)
        solver.beta_lin = one(T)
    end

    @inbounds if predict
        # solution update
        for i = 1:solver.model.nlmi
            solver.Xn[i] .= solver.X[i] .+ solver.alpha[i] .* solver.delX[i]
            solver.Sn[i] .= solver.S[i] .+ solver.beta[i] .* solver.delS[i]
            RNTᵢ, Dᵢ = solver.RNT[i], solver.D[i]
            mul!(RNTᵢ, solver.delS[i], solver.G[i])
            tmp = mul!(solver.delS[i], solver.delX[i], RNTᵢ)
            mul!(RNTᵢ, solver.Gi[i], tmp)
            for j in 1:size(RNTᵢ, 2), k in 1:j
                RNTᵢ[j, k] = RNTᵢ[k, j] = -(RNTᵢ[k, j] + RNTᵢ[j, k]) / (Dᵢ[j] + Dᵢ[k])
            end
        end
    else
        αₘ = min(minimum(solver.alpha, init=typemax(T)), solver.alpha_lin)
        βₘ = min(minimum(solver.beta, init=typemax(T)), solver.beta_lin)
        solver.y .+= βₘ .* solver.dely
        for i in 1:solver.model.nlmi
            axpy!(αₘ, solver.delX[i], solver.X[i])
            axpy!(βₘ, solver.delS[i], solver.S[i])
        end
    end

    return
end

function find_step_lin!(solver::Solver{T}, predict::Bool) where {T}
    copyto!(solver.delS_lin, solver.Rd_lin)
    mul!(solver.delS_lin, solver.model.A_lin', solver.dely, -one(T), true)
    if predict
        solver.delX_lin .= solver.X_lin .* (-one(T) .- solver.Si_lin .* solver.delS_lin)
    else
        solver.delX_lin .= solver.X_lin .* (-one(T) .- solver.Si_lin .* solver.delS_lin) .+
                           (solver.sigma * solver.mu) .* solver.Si_lin .+ solver.RNT_lin
    end
    mimiX_lin = minimum(splat(/), zip(solver.delX_lin, solver.X_lin), init=typemax(T))
    if mimiX_lin > T(-1//1_000_000)
        solver.alpha_lin = T(99//100)
    else
        solver.alpha_lin = min(one(T), -solver.tau / mimiX_lin)
    end
    mimiS_lin = minimum(splat(/), zip(solver.delS_lin, solver.S_lin), init=typemax(T))
    if mimiS_lin > T(-1//1_000_000)
        solver.beta_lin = T(99//100)
    else
        solver.beta_lin = min(one(T), -solver.tau / mimiS_lin)
    end

    if predict
        # solution update
        solver.Xn_lin .= solver.X_lin .+ solver.alpha_lin .* solver.delX_lin
        solver.Sn_lin .= solver.S_lin .+ solver.beta_lin .* solver.delS_lin

        solver.RNT_lin .= .-(solver.delX_lin .* solver.delS_lin) .* solver.Si_lin
    else
        solver.X_lin .+= min(minimum(solver.alpha, init=typemax(T)), solver.alpha_lin) .* solver.delX_lin
        solver.S_lin .+= min(minimum(solver.beta, init=typemax(T)), solver.beta_lin) .* solver.delS_lin
        solver.S_lin_inv .= inv.(solver.S_lin)
    end

    return
end
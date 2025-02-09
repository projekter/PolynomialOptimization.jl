struct PreconditionerNone end

(t::PreconditionerNone)(Mx::Vector{T}, x::Vector{T}) where {T} = copyto!(Mx, x)

prepare_prec!(solver::Solver, ::PreconditionerNone) = nothing

mutable struct PreconditionerAlpha{T<:Real,I<:Integer,AT<:Union{UniformScaling{T},EfficientCholmod{T,CHOLMOD.Factor{T,I}}}}
    const solver::Solver{T,I}
    const ntot::Int
    const Umat::Vector{Matrix{T}}
    const Z::Vector{LowerTriangular{T,Matrix{T}}}
    const cholS::Cholesky{T,Matrix{T}}
    const y22::Vector{Vector{T}}
    const y33::Vector{T}
    const y33s::Vector{Base.ReshapedArray{T,2,SubArray{T,1,Vector{T},Tuple{UnitRange{Int}},true},Tuple{}}}
    const AAinvx::Vector{T}
    const t::Matrix{T}
    AATtau::AT

    function PreconditionerAlpha(solver::Solver{T,I}) where {T<:Real,I<:Integer}
        ntot = sum(solver.model.coneDims, init=0)
        entot = solver.erank * ntot
        y33 = Vector{T}(undef, entot)
        y33s = Vector{Base.ReshapedArray{T,2,SubArray{T,1,Vector{T},Tuple{UnitRange{Int}},true},Tuple{}}}(
            undef, solver.model.nlmi
        )
        ii = 1
        @inbounds for (ilmi, n) in enumerate(solver.model.coneDims)
            dim = n * solver.erank
            y33s[ilmi] = reshape(@view(y33[ii:ii+dim-1]), n, solver.erank)
            ii += dim
        end
        new{T,I,iszero(solver.model.nlin) ? UniformScaling{T} : EfficientCholmod{T,CHOLMOD.Factor{T,I}}}(
            solver,
            ntot,
            [Matrix{T}(undef, s, solver.erank) for s in solver.model.coneDims],
            [LowerTriangular(Matrix{T}(undef, s, s)) for s in solver.model.coneDims],
            Cholesky(Matrix{T}(undef, entot, entot), 'U', 0),
            [Vector{T}(undef, s^2) for s in solver.model.coneDims],
            y33, y33s,
            Vector{T}(undef, solver.model.n),
            Matrix{T}(undef, solver.model.n, solver.erank * ntot)
        )
    end
end

function (t::PreconditionerAlpha{T})(Mx::Vector{T}, x::Vector{T}) where {T<:Real}
    nvar = length(x)
    nlmi = length(t.solver.model.A)

    AAinvx = ldiv!(t.AATtau, copyto!(t.AAinvx, x))

    Threads.@threads for ilmi in 1:nlmi
        @inbounds begin
            n = size(t.Umat[ilmi], 1)
            yy = t.y22[ilmi]
            yym = reshape(yy, n, n)
            mul!(yy, t.solver.model.A[ilmi]', AAinvx)
            lmul!(t.Z[ilmi]', yym)
            mul!(t.y33s[ilmi], yym, t.Umat[ilmi])
        end
    end
    ldiv!(t.cholS, t.y33)

    fill!(Mx, zero(T))
    Threads.@threads for ilmi in 1:nlmi
        @inbounds begin
            n = size(t.Umat[ilmi], 1)
            yym = reshape(fill!(t.y22[ilmi], zero(T)), n, n)
            for (Ucol, target) in zip(eachcol(t.Umat[ilmi]), eachcol(t.y33s[ilmi]))
                lmul!(t.Z[ilmi], target)
                BLAS.ger!(one(T), target, Ucol, yym)
            end
        end
    end
    @inbounds for ilmi in 1:nlmi
        mul!(Mx, t.solver.model.A[ilmi], t.y22[ilmi], true, true)
    end
    ldiv!(t.AATtau, Mx)

    axpby!(true, AAinvx, -one(T), Mx)
    return Mx
end

function prepare_prec!(solver::Solver{T}, preconditioner::PreconditionerAlpha{T,<:Integer,AT}) where {T<:Real,AT}
    nlmi = solver.model.nlmi
    iszero(nlmi) && return
    k = solver.erank
    nvar = solver.model.n
    ntot = preconditioner.ntot

    lbt = 1
    δ = zero(T)
    @inbounds for ilmi in 1:nlmi
        lambdaf, vectf = eigen(Symmetric(solver.W[ilmi], :U))
        vect_l = @view(vectf[:, end-k+1:end])
        lambda_l = @view(lambdaf[end-k+1:end])
        vect_s = @view(vectf[:, 1:end-k])
        lambda_s = @view(lambdaf[1:end-k])

        ttau = lambda_s[1]
        if solver.aamat != AMAT_AᵀA
            ttau = (ttau + sum(lambda_s, init=zero(T)) / length(lambda_s)) / 2 - T(1//10_000_000_000_000)
        end
        lambda_l .= sqrt.(lambda_l .- ttau)
        mul!(preconditioner.Umat[ilmi], vect_l, Diagonal(lambda_l))

        fill!(lambda_l, ttau)
        lambdaf .= sqrt.(max.(zero(T), lambdaf)) # W is supposed to be positive definite, so the max is just for numerics
        rmul!(vectf, Diagonal(lambdaf))
        W0 = BLAS.syrk!('L', 'N', true, vectf, false, parent(preconditioner.Z[ilmi]))
        BLAS.syrk!('L', 'N', true, preconditioner.Umat[ilmi], T(2), W0)
        cholesky!(Symmetric(W0, :L), check=false)

        δ += ttau^2
    end

    if AT <: UniformScaling{T}
        preconditioner.AATtau = δ * I
    else
        # TODO: make it efficient
        AATtau = solver.model.A_lin * Diagonal(solver.X_lin .* solver.S_lin_inv) * solver.model.A_lin' + δ * I
    end

    t = preconditioner.t
    @inbounds if k > 1 #slow formula
        for ilmi in 1:nlmi
            n = size(solver.W[ilmi], 1)
            TT = kron(preconditioner.Umat[ilmi], preconditioner.Z[ilmi]) # TODO: lazy kron
            mul!(@view(t[:, lbt:lbt+k*n-1]), solver.model.A[ilmi], TT)
            lbt += k * n
        end
        if AT <: UniformScaling{T}
            lmul!(sqrt(preconditioner.AATtau), t)
        else
            preconditioner.AATtau = EfficientCholmod(cholesky(Symmetric(AATtau, :U), check=false))
            lmul!(preconditioner.AATtau.U, t)
        end
    else # fast formula
        if !(AT <: UniformScaling{T})
            AATtau_d = sqrt.(inv.(diag(AATtau)))
            preconditioner.AATtau = EfficientCholmod(cholesky(Symmetric(AATtau, :U), check=false))
        end
        for ilmi in 1:solver.model.nlmi
            n = size(solver.W[ilmi], 1)

            ii_, jj_, aa_ = findnz(solver.model.A[ilmi])
            jj_ .-= one(eltype(jj_))
            qq_ = similar(jj_)
            pp_ = jj_
            Umatᵢ = preconditioner.Umat[ilmi]
            for i in eachindex(jj_)
                qq_[i], pp_[i] = divrem(jj_[i], n) .+ 1
                if AT <: UniformScaling{T}
                    aa_[i] *= Umatᵢ[qq_[i]]
                else
                    aa_[i] *= Umatᵢ[qq_[i]] * AATtau_d[ii_[i]]
                end
            end
            if AT <: UniformScaling{T}
                rmul!(aa_, inv(sqrt(preconditioner.AATtau[1, 1])))
            end
            AU = SparseArrays.sparse!(ii_, pp_, aa_, nvar, n)
            mul!(@view(t[1:nvar, lbt:lbt+k*n-1]), AU, preconditioner.Z[ilmi])
            lbt = lbt + k * n
        end
    end
    S = BLAS.syrk!('U', 'T', true, t, false, preconditioner.cholS.factors)
    # Schur complement for the SMW formula
    @inbounds for i in diagind(S)
        S[i] += one(T)
    end
    cholesky!(Symmetric(S, :U), check=false)
    return
end

mutable struct PreconditionerBeta{T<:Real,I<:Integer,AT<:Union{T,Vector{T}}}
    const solver::Solver{T,I}
    const ntot::Int
    AATtau::AT

    PreconditionerBeta(solver::Solver{T,I}) where {T<:Real,I<:Integer} =
        new{T,I,iszero(solver.model.nlin) ? T : Vector{T}}(
            solver,
            sum(solver.model.coneDims, init=0),
            iszero(solver.model.nlin) ? zero(T) : Vector{T}(undef, solver.model.n)
        )
end

function (t::PreconditionerBeta{T})(Mx::Vector{T}, x::Vector{T}) where {T<:Real}
    copyto!(Mx, x)
    Mx ./= t.AATtau
    return Mx
end

function prepare_prec!(solver::Solver{T}, preconditioner::PreconditionerBeta{T,<:Integer,AT}) where {T<:Real,AT}
    nlmi = solver.model.nlmi
    iszero(nlmi) && return
    k = solver.erank
    nvar = solver.model.n
    ntot = preconditioner.ntot

    δ = zero(T)
    @inbounds for ilmi in 1:nlmi
        lambda_s = eigvals(Symmetric(solver.W[ilmi], :U), 1:size(solver.W[ilmi], 1)-k)

        ttau = lambda_s[1]
        if solver.aamat != AMAT_AᵀA
            ttau = (ttau + sum(lambda_s, init=zero(T)) / length(lambda_s)) / 2 - T(1//10_000_000_000_000)
        end
        δ += ttau^2
    end
    if AT <: UniformScaling
        preconditioner.AATtau = δ
    else
        fill!(preconditioner.AATtau, δ)
        @inbounds for i in 1:solver.model.nlin
            α = solver.X_lin[i] * solver.S_lin_inv[i]
            v = @view(solver.model.A_lin[:, i])
            for (j, v) in zip(rowvals(v), nonzeros(v))
                preconditioner.AATtau[j] += α * v^2
            end
        end
    end

    return
end

# union type instead of abstract to allow for faster dispatch
const PreconditionerUnion{T<:Real,I<:Integer} = Union{PreconditionerNone,<:PreconditionerAlpha{T,I},<:PreconditionerBeta{T,I}}
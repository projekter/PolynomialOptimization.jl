# This is an implementation of the SketchyCGAL solver, https://doi.org/10.1137/19M1305045, extended to work with multiple
# semidefinite matrices of varying size. This cannot simply be reproduced by using block-diagonals and the usual algorithm,
# since the block-diagonal constraints actually increase the rank of the full matrix, which defeats the purpose of a low-rank
# solver.
module SketchyCGAL

using IterativeSolvers, LinearAlgebra, LinearMaps, PositiveFactorizations, Printf, Random, SparseArrays
using ...PolynomialOptimization: @assert, @inbounds, @verbose_info

export Status, sketchy_cgal

"""
    Status{R}

This struct contains the current information for the Sketchy CGAL solver. Per-iteration callbacks will receive this structure
to gather current information.

See also [`sketchy_cgal`](@ref).
"""
mutable struct Status{R}
    # static information
    max_iter::UInt
    time_limit::UInt
    ϵ::R
    # dynamic information
    iteration::UInt
    time::UInt
    objective::R
    infeasibility::R
    infeasibility_rel::R
    suboptimality::R
    suboptimality_rel::R

    Status(max_iter::Integer, time_limit::Integer, ϵ::R) where {R} = new{R}(UInt(max_iter), UInt(time_limit), ϵ)
end

@doc raw"""
    sketchy_cgal(A, b, C; α=(0, 1), rank, ϵ=1e-4, max_iter=0, time_limit=0, verbose=false, β₀=1, K=∞, method=:auto,
        callback=(_) -> ()[, A_norm][, A_normsquare])

Enhanced implementation of the [SketchyCGAL algorithm](https://doi.org/10.1137/19M1305045). This solves the following problem:
```math
    \min \bigl\{
        \sum_i \operatorname{tr}(C_i X_i) :
            \sum_i \operatorname{tr}(A_{j, i} X_i) = b_j \ \forall j,
            \alpha_1 \leq \sum_i \operatorname{tr}(X_i) \leq \alpha_2,
            X_i \succeq 0 \ \forall i
    \bigr\}
```
The function returns the optimization status, objective value, and the optimal matrix ``X`` (in the form of an `Eigen`
factorization object).

## Parameters
- `A` is an `AbstractMatrix` whose entries are of type `AbstractMatrix` themselves. Alternatively, `A` can also be an
  `AbstractVector` of `AbstractMatrices`; in this case, ``A_{j, i}`` is given by taking the matrix `A[i]` and reshaping its
  columns into a square matrix.
- `b` is an `AbstractVector` of real numbers
- `C` is an `AbstractVector` whose entries are of type `AbstractMatrix` themselves
- `α` is a 2-tuples of nonnegative numbers, where the numbers defines the bounds on the sum of all traces
- `rank` controls the rank that is used to approximate the end result. It might either be an integer, which then puts the same
  rank constraint on all ``X_i``, or a tuple/vector of integers, which allows to specify different rank constraints.
- The solution accuracy can be controlled by the parameter `ϵ`; however, no more than `max_iter` iterations are carried out,
  and no more iterations will be performed if `time_limit` was exceeded (in seconds), regardless of `ϵ`. Set any of those three
  parameters to zero to disable the check.
- The `callback` may be called after each iteration and will receive a [`Status`](@ref) as parameter. If the callback returns
  `false`, the iteration will be the last one.
- The parameters `β₀` and `K` allow to tune the optimization. `β₀` is a smoothing, `K` limits the dual vector to a generalized
  sphere of radius `K` around the origin.
- The `method` determines the way in which the smallest eigenvalue and its eigenvector are calculated during each iteration.
  Possible values are (this might also be a vector/tuple, to specify the method for each ``X_i``)
  - `:lanczos_space` (uses the space-efficient implementation described in the SketchyCGAL paper, memory is linear in the
    problem size)
  - `:lanczos_time` (uses the same principle, but can save about half of the operations by using more memory: quadratic in the
    problem size)
  - `:lobpcg_fast` (uses the LOBPCG solver from the [`IterativeSolvers`](https://iterativesolvers.julialinearalgebra.org/dev/)
    package, bounding the number of iterations with the same heuristic as for the Lanczos methods)
  - `:lobpcg_accurate` (bounds the error instead to `ϵ/100`)
  - `:auto` chooses `:lobpcg_accurate` for problem sizes smaller than 10, `:lanczos_time` for problem sizes less than 11500
    (where roughly 1 GiB is required for the speedup), and `:lanczos_space` for all larger problems.
- `A_normsquare` (or `A_norm`) is supposed to hold the sum of the squares of the Frobenius-to-ℓ₂-norms of all the linear
  operators contained in the columns of `A`. If both parameters are omitted, it is calculated automatically; however, this
  requires memory that scales quartically in the largest side dimension of the `A` (and may not be supported for all
  `AbstractMatrix` types). If both parameters are specified and their values are not consistent, the behavior is undefined.

## Optimization status values
- `:optimal`: the desired accuracy `ϵ` was reached in both the relative suboptimality gap as well as the relative infeasibility
- `:max_iter`: the maximum number of iterations was reached
- `:timeout`: the maximum computation time was hit
- `:canceled`: the callback returned `false`
- `:unknown`: an internal error has happened

See also [`Status`](@ref).



    sketchy_cgal(primitive1!, primitive2!, primitive3!, n, b, α; rank, primitive3_norm=0, primitive3_normsquare=0, ϵ=1e-8,
        max_iter=10_000, verbose=false, rescale_C=1, rescale_A=[1, ...], rescale_X=1, β₀=1, K=∞, method=:auto,
        callback=(_) -> ())

This is the black-box version that allows for matrix-free operations. `n` is a tuple or vector that indicates the side
dimensions of the semidefinite variables, `b` is the right-hand side of the constraint vector (of length `d`); and the
primitives effectively calculate
- `primitive1!(v, u, i, α, β) = (v = α * C[i] * u + β * v)`, ``u, v \in \mathbb R^{n_i}``
- `primitive2!(v, u, z, i, α, β) = (v = α * adjoint(A[:, i])(z) * u + β * v)`, ``u, v \in \mathbb R^{n_i}``,
  ``z \in \mathbb R^d``
- `primitive3!(v, u, i, α, β) = (v = α * A[:, i](u * u') + β * v)`, ``u \in \mathbb R^{n_i}``, ``v \in \mathbb R^d``
`A[:, i](X)` is the linear map `[⟨X, A[1, i]⟩, ..., ⟨X, A[d, i]⟩]` and `adjoint(A[:, i])(z) = sum(z[j] * A[j, i])`. All of them
must also return their outputs (which is `v`).

If you are able to calculate these oracles faster or more memory-efficiently than the straightforward implementation (which is
based on `mul!`), use the blackbox method.
It is recommended to obey the following normalization conditions:
```math
    \sum_i lVert C_i\rVert_{\mathrm F}^2 = 1;
    \quad
    \sum_i \lVert primitive3!_i\rVert_{\mathrm F \to \ell_2}^2 = 1;
    \quad
    \sum_i \lVert A_{1, i}\rVert^2 = \sum_i \lVert A_{2, i}\rVert^2 = \dotsb
```
However, in any case, you need to specify the norm of `primitive3!` (i.e., the supremum of `norm(primitive3!)` applied to
matrices with unit Frobenius norm) in the parameter `primitive3_norm` (or `primitive3_normsquare`, which is the sum of the
individual normsquares). If the norm is unavailable, you need to at least give a lower bound. You may not specify both
parameters inconsistently, else the behavior is undefined.
If you had to rescale those matrices in order to achieve this normalization condition, you may pass the corresponding rescaling
factors in `rescale_C` (implying that all `C` have been scaled by this one factor) and `rescale_A` (implying that a whole row
of `A` has been scaled by one element from this vector).
Additionally, the upper bound in `α` should be one for a better performance, which is achievable through `rescale_X` (implying
that all `X` have been scaled by this factor). These are posterior factors that indicate the multiplication that has been done
before calling the function in order to enforce compliance. They will be taken into account when calculating the termination
criteria (such that `ϵ` then corresponds to the original problem and not the rescaled one), filling the status structure or
verbose output, and the final values of primal objective and optimal ``X``.
"""
function sketchy_cgal(primitive1!, primitive2!, primitive3!,
    @nospecialize(n::Union{Integer,Tuple{Integer},AbstractVector{<:Integer}}), b::AbstractVector{R}; α::Tuple{R,R},
    @nospecialize(rank::Union{Integer,Tuple{Integer},AbstractVector{<:Integer}}), primitive3_norm::R=zero(R),
    primitive3_normsquare::R=zero(R), ϵ::R=1e-4, max_iter::Integer=10_000, time_limit::Integer=0, rescale_C::R=one(R),
    rescale_A::Vector{R}=fill(one(R), length(b)), rescale_X::R=one(R), β₀::R=one(R), K::R=R(Inf), verbose::Bool=false,
    callback::Union{Nothing,Function}=nothing,
    @nospecialize(method::Union{Symbol,Tuple{Symbol},AbstractVector{<:Symbol}}=:auto)) where {R<:Real}
    (ϵ ≥ zero(R) && max_iter ≥ 0 && time_limit ≥ 0 && (primitive3_norm > zero(R) || primitive3_normsquare > zero(R)) &&
        α[1] ≥ zero(R) && α[2] ≥ α[1] && α[2] > zero(R)) ||
        throw(ArgumentError("Tolerances and limits must be nonnegative and well-ordered."))
    (ϵ > zero(R) || max_iter > 0 || time_limit > 0) ||
        throw(ArgumentError("At least one of the following must be strictly positive: tolerance, iteration limit, time limit."))
    if n isa Integer
        n = [n]
    end
    N = length(n)
    if rank isa Integer
        rank ≥ 1 || throw(ArgumentError("The rank must exceed 1."))
        rank = fill(rank, N)
    else
        (length(rank) == N && all(rank .≥ 1)) ||
            throw(ArgumentError("Each semidefinite variable must have a rank at least 1."))
    end
    if method isa Symbol
        methods = fill(method, N)
    else
        length(method) == N || throw(ArgumentError("Each semidefinite variable must have a method."))
        methods = method isa Tuple ? [method...] :
                                     (method isa Vector{Symbol} ? method : collect(method)::Vector{Symbol})
    end
    @inbounds for i in eachindex(methods)
        if methods[i] === :auto
            if n[i] < 10
                methods[i] = :lobpcg_accurate
            elseif n[i] < 11_500
                methods[i] = :lanczos_time
            else
                methods[i] = :lanczos_space
            end
        else
            methods[i] ∈ (:lanczos_space, :lanczos_time, :lobpcg_fast, :lobpcg_accurate) ||
                throw(ArgumentError("Method $i had an invalid value. Possible methods are `:lanczos_space`, `:lanczos_time`, \
                                     `:lobpcg_fast`, `: lobpcg_accurate`."))
        end
    end
    starting_time = time_ns()
    info = Status(max_iter, time_limit, ϵ)
    detail_calc = verbose || !isnothing(callback)
    if !iszero(time_limit)
        time_limit = starting_time + 1_000_000_000 * time_limit
    end
    @verbose_info("SketchyCGAL solver, implemented by PolynomialOptimization, using eigensolver method(s) ", union(methods),
        "\n", "Iteration | Primal objective | Suboptimality | Infeasibility | rel. subopt. | rel. infeas. | Time")

    d = length(b)
    sqrtα = sqrt.(α)
    status = :unknown
    # Scale problem data (or better: problem data is scaled, but adjust termination criteria appropriately)
    obj_rescale = 1 / (rescale_C * rescale_X)
    infeas_rescale = (1 / rescale_X) ./ rescale_A
    infeas_rel_factor = 1 / max(sqrt(sum(x -> prod(x, init=one(R)::R)^2, zip(b, infeas_rescale), init=zero(R))::R), one(R))
    if primitive3_norm > 0
        primitive3_normsquare = primitive3_norm^2
    end
    # β₀ ← 1 and K ← ∞ (we expose these parameters)
    # NystromSketch.Init(n, R)
    sketch = NystromSketch{R}(n, rank)
    # z ← 0_d and y ← 0_d
    z = zeros(R, d, N)
    y = zeros(R, d)
    # we also define some temporaries to avoid allocations and speed up calculations in the loop
    ∑zminusb = -b
    tmpd = Vector{R}(undef, d)
    tmpv = let data = Vector{R}(undef, maximum(n))
        @inbounds [@view(data[1:nᵢ]) for nᵢ in n]
    end
    eig_tmp = setup_approx_min_evec.(Val.(methods), n, ϵ, (tmpd,), primitive1!, primitive2!)
    ξ = Vector{R}(undef, N)
    v = Vector{Vector{R}}(undef, N)
    p = zero(R)
    trace = zero(R)
    # for t ← 1, 2, 3, ..., T do
    t = 1
    @inbounds while true
        # β ← β₀ * √(t+1) and η ← 2/(t+1)
        β = β₀ * sqrt(t +1)
        η = 2 / (t +1)
        # [ξ, v] ← ApproxMinEvec(C + A*(y + β(z - b)); q) with q = t^(1/4) log n
        tmpd .= y .+ β .* ∑zminusb
        for i in 1:N
            ξ[i], v[i] = approx_min_evec(t, i, eig_tmp[i])::Tuple{R,Vector{R}}
        end
        # we are now interested in the smallest eigenvalue of all
        imin = argmin(ξ)
        ξmin, vmin = ξ[imin], v[imin]
        # depending on the sign of ξ, we need to pick either the upper or lower trace bound.
        # Note that here, we already multiply v by sqrtα. We need to do this anyway in the end for the sketch update, so we can
        # just do it here and avoid some α multiplications afterwards.
        lmul!(ξmin > 0 ? sqrtα[1] : sqrtα[2], vmin)
        # check the stopping criterion (for the current point, where we need ξ)
        finish = false
        if !iszero(ϵ) || detail_calc
            # the infeasibility is just z - b
            info.infeasibility = sqrt(sum(x -> prod(x, init=one(R))^2, zip(∑zminusb, infeas_rescale), init=zero(R))::R)
            info.infeasibility_rel = info.infeasibility * infeas_rel_factor
            stop_feasible = info.infeasibility_rel ≤ ϵ
            if stop_feasible || detail_calc
                # we can bound the suboptimality via
                # g - ⟨y, z - b⟩ - 1/2 β ‖ z - b ‖^2
                # where g = p + ⟨y + β (z - b), z⟩ - λₘᵢₙ(D)
                # and λₘᵢₙ(D) ≈ ξ
                info.suboptimality = (p + sum(dot(tmpd, zᵢ) for zᵢ in _eachcol(z), init=zero(R)) - ξmin - dot(y, ∑zminusb) -
                    β * norm(∑zminusb)^2/2) * obj_rescale
                info.suboptimality_rel = info.suboptimality / max(abs(p * obj_rescale), 1)
                finish = stop_feasible && info.suboptimality_rel ≤ ϵ
            end
        end
        if detail_calc
            info.iteration = t
            info.objective = p * obj_rescale
            info.time = (time_ns() - starting_time) ÷ 1_000_000_000
            if verbose && (finish || isone(t) || iszero(t % 100))
                @verbose_info(@sprintf("%9d | %16g | %13g | %13g | %12g | %12g | %02d:%02d",
                    t, # iteration
                    info.objective, # primal objective
                    info.suboptimality, # actual upper bound to suboptimality
                    info.infeasibility, # actual infeasibility
                    info.suboptimality_rel, # relative value of suboptimality, relevant for stopping
                    info.infeasibility_rel, # relative value of infeasibility, relevant for stopping
                    info.time ÷ 60, info.time % 60))
            end
            if !isnothing(callback) && callback(info) === false
                status = :canceled
                break
            end
        end
        if finish
            status = :optimal
            break
        end
        if !iszero(max_iter) && t == max_iter
            status = :max_iter
            break
        end
        if !iszero(time_limit) && time_ns() ≥ time_limit
            status = :timeout
            break
        end

        # update pₜ₊₁ = (1 - η) pₜ + η * α * vₜ*(C vₜ) - strictly speaking, we don't need this for zero ϵ, but we want to keep
        # track of the objective anyway (and even for no output, we output it in the end)
        p = (1 - η) * p + η * dot(vmin, primitive1!(tmpv[imin], vmin, imin, true, false))
        # z = (1 - η) z + η A(α v v*)
        primitive3!(@view(z[:, imin]), vmin, imin, η, one(R) - η)
        sum!(∑zminusb, z)
        ∑zminusb .-= b
        trace = (1 - η) * trace + η * (ξmin > 0 ? α[1] : α[2])
        # γ is the largest solution to γ ‖ zₜ₊₁ - b ‖² ≤ 4α²β₀ / (t +1)^(3/2) ‖ A ‖^2 and 0 ≤ γ ≤ β₀,
        # which can also be written as γ ‖ zₜ₊₁ - b ‖² ≤ β η²α² ‖ A ‖² and 0 ≤ γ ≤ β₀
        # If needed, set γ = 0 to prevent ‖ yₜ₊₁ ‖ > K (which "in practice is not necessary").
        γ = min(β * (η * α[2])^2 * primitive3_normsquare / norm(∑zminusb)^2, β₀)
        # Note that in MATLAB, this strangely uses β₊ = β₀ * sqrt(t +2) (so the β for the next iteration), but still η from
        # this iteration...
        # yₜ₊₁ = yₜ + γ (zₜ₊₁ - b)
        tmpd .= y .+ γ .* ∑zminusb
        if norm(tmpd) ≤ K
            y .= tmpd
        end
        # NystromSketch.RankOneUpdate(√α v, η)
        rank_one_update!(sketch, vmin, imin, η)

        t += 1
    end
    # [U, Λ] ← NystromSketch.Reconstruct()
    U, Λ = reconstruct(sketch)
    # Λ ← Λ + (α - tr(Λ)) I / R
    # Here, we do not use α as the trace reference, since α might actually define an interval. Instead, we kept track of what
    # the trace was supposed to be (assuming an ideal storage). We then compute the actual trace based on the sketch
    # reconstruction and perform the correction step based on this.
    trace_correction = (trace - sum(sum, Λ, init=zero(R))::R) / sum(rank, init=zero(R))::R
    # ^ This way of correction corresponds to the correction that would be done if all of the matrices had been assembled in a
    # large block matrix.
    for Λᵢ in Λ # rank == length(Λᵢ)
        Λᵢ .= (Λᵢ .+ trace_correction) ./ rescale_X
    end
    # X = U Λ U*. Instead of returning the full matrix that arises in this way, we give back an appropriate factorization.
    return status, p * obj_rescale, Eigen.(Λ, U)
end

struct OpNorm{V}
    v::V
end

function (o::OpNorm)(y, x)
    # m = ∑ᵢ o.v[i] o.v[i]ᵀ
    # → m*x = ∑ᵢ dot(o.v[i], x) o.v[i]
    firstel, rest = Iterators.peel(o.v)
    y .= dot(firstel, x) .* firstel
    for vᵢ in rest
        axpy!(dot(vᵢ, x), vᵢ, y)
    end
    return y
end

_eachrow(x) = eachrow(x)
_eachrow(x::Transpose) = _eachcol(parent(x))
_eachcol(x) = eachcol(x)
_eachcol(x::Transpose) = _eachrow(parent(x))
Base.@propagate_inbounds _view(x, a, b) = view(x, a, b)
Base.@propagate_inbounds _view(x::Transpose, a, b) = view(parent(x), b, a)

@inline function _pseudodot(A::AbstractVector{T}, y::AbstractVector{T}) where {T<:Real}
    # calculates dot(y, reshape(A)[:, :], y) where the single dimension of the matrix A is split in two
    @boundscheck (length(A) == length(y)^2) || throw(DimensionMismatch())
    s = zero(T)
    y₁ = firstindex(y)
    ia = firstindex(A)
    @inbounds for y₁ in y
        @simd for y₂ in eachindex(y)
            s += A[ia+y₂-firstindex(y₂)] * y₁ * y[y₂]
        end
        ia += length(y)
    end
    return s
end

@inline function _pseudodot(A::SparseArrays.SparseVectorUnion{T}, y::AbstractVector{T}) where {T<:Real}
    @boundscheck (length(A) == length(y)^2) || throw(DimensionMismatch())
    s = zero(T)
    dim = length(y)
    y₁ = firstindex(y)
    lastincol = dim
    for (i, v) in zip(rowvals(A), nonzeros(A))
        while i > lastincol # avoid division, assuming the values are not too widely spaced
            y₁ += 1
            lastincol += dim
        end
        y₂ = i - (lastincol - dim +1) + 1 + firstindex(y) - firstindex(axes(A, 1))
        s += v * y[y₁] * y[y₂]
    end
    return s
end

function sketchy_cgal(A::Union{<:AbstractMatrix{<:AbstractMatrix{R}},<:AbstractVector{<:AbstractMatrix{R}}},
    b::AbstractVector{R}, C::AbstractVector{<:AbstractMatrix{R}}; verbose::Bool=false, kwargs...) where {R<:Real}
    n = LinearAlgebra.checksquare.(C)
    N = length(n)
    if A isa AbstractMatrix
        (size(A, 1) == length(b) && size(A, 2) == length(C) == N && size(A, 1) ≥ 1) ||
            throw(ArgumentError("Invalid matrix dimensions."))
        for (nᵢ, Acol) in zip(n, _eachcol(A))
            all(x -> LinearAlgebra.checksquare(x[1]) == nᵢ, Acol) || throw(ArgumentError("Invalid matrix dimensions."))
        end
    else
        length(A) == length(C) == N || throw(ArgumentError("Invalid matrix dimensions."))
        for (nᵢ, Aᵢ) in zip(n, A)
            size(Aᵢ) == (length(b), nᵢ^2) || throw(ArgumentError("Invalid matrix dimensions."))
        end
    end
    kwargs = Dict(kwargs) # kwargs are immutable

    d = length(b)
    @verbose_info("Calculating the rescaling parameters")
    rescale_C = 1 / sqrt(sum(LinearAlgebra.norm_sqr, C, init=zero(R)))
    # now we rescale the A such that their Frobenius norms are all the same. Note that here, we need to combine the As of one
    # column together
    rescale_A = Vector{R}(undef, d)
    @inbounds rescale_A[1] = one(R)
    if A isa AbstractMatrix
        let target=sqrt(sum(LinearAlgebra.norm_sqr, first(_eachrow(A)), init=zero(R))::R)
            for (j, Arow) in Iterators.drop(enumerate(_eachrow(A)), 1)
                @inbounds rescale_A[j] = target / sqrt(sum(LinearAlgebra.norm_sqr, Arow, init=zero(R))::R)
            end
        end
    else
        # TODO: if we have SparseMatrixCSC, this can be optimized by traveling all rowvals and incrementing rescale_A for this
        # row, then in the end inverting.
        let target=sqrt(sum(∘(LinearAlgebra.norm_sqr, first, _eachrow), A, init=zero(R))::R)
            for j in 2:d
                @inbounds rescale_A[j] = target / sqrt(sum(x -> LinearAlgebra.norm_sqr(_view(x, j, :)), A, init=zero(R))::R)
            end
        end
    end
    # However, we must also take care of the operator norms of all the A, which must be one. The norm calculation is a bit
    # costly, but at least it's doable. The user may want to provide some insight and already prespecify the norm.
    if !haskey(kwargs, :A_normsquare) && !haskey(kwargs, :A_norm)
        @verbose_info("Determining operator norm of A")
        A_normsquare = zero(R)
        powϵ = convert(R, get(kwargs, :ϵ, zero(R)))::R
        if iszero(powϵ)
            powϵ = R(1//1000)
        end
        @inbounds for j in 1:N
            A_normsquare += real(
                powm(LinearMap{R}(
                    OpNorm(A isa AbstractMatrix ? map(vec, _view(A, :, j)) : _eachrow(A[j])),
                    n[j]^2, issymmetric=true, isposdef=true, ismutating=true), tol=powϵ, verbose=n[j] > 500 && verbose
                )[1]
            )
        end
        @verbose_info("Found operator norm: ", sqrt(A_normsquare))
        # And after we got the norm, we have to rescale everything once again.
        rescale_A ./= sqrt(A_normsquare)
    elseif haskey(kwargs, :A_normsquare)
        A_normsquare = convert(R, pop!(kwargs, :A_normsquare))::R
    else
        A_normsquare = (convert(R, pop!(kwargs, :A_norm))::R)^2
    end
    #@verbose_info("Rescale C by $rescale_C, A by $rescale_A.") # may be pretty long

    alpha = pop!(kwargs, :α, (zero(R), one(R)))::Tuple{R,R}
    @inbounds rescale_X = 1 / alpha[2]
    status, obj, X = @inbounds sketchy_cgal(
        (v, u, i, α, β) -> mul!(v, C[i], u, α * rescale_C, β),
        if A isa AbstractMatrix
            (v, u, z, i, α, β) -> begin
                mul!(v, A[1, i], u, α * z[1] * rescale_A[1], β)
                for j in 2:d
                    mul!(v, A[j, i], u, α * z[j] * rescale_A[j], one(R))
                end
                v
            end
        else
            let tmp=[let q=Matrix{R}(undef, nᵢ, nᵢ); (q, vec(q)) end for nᵢ in n]
                (v, u, z, i, α, β) -> begin
                    if isone(i)
                        z .*= rescale_A
                    end
                    mul!(tmp[i][2], transpose(A[i]), z)
                    mul!(v, tmp[i][1], u, α, β)
                    if i == N
                        z ./= rescale_A
                    end
                    v
                end
            end
        end,
        if A isa AbstractMatrix
            (v, u, i, α, β) -> begin
                for j in 1:d
                    v[j] = α * rescale_A[j] * dot(u, A[j, i], u) + β * v[j]
                end
                v
            end
        else
            (v, u, i, α, β) -> begin
                for (j, Aⱼ) in enumerate(_eachrow(A[i]))
                    v[j] = α * rescale_A[j] * _pseudodot(Aⱼ, u) + β * v[j]
                end
                v
            end
        end,
        n, b .* rescale_A .* rescale_X; α=alpha .* rescale_X,
        verbose, rescale_C, rescale_A, rescale_X, primitive3_normsquare=one(R),
        kwargs...
    )
    return status, obj, X
end

struct LanczosTimeTmp{R,P1,P2}
    n::Int
    logn::R
    vmat::Matrix{R}
    ρ::Vector{R}
    ω::Vector{R}
    v::Vector{R}
    tmpd::Vector{R}
    primitive1!::P1
    primitive2!::P2
end

function setup_approx_min_evec(::Val{:lanczos_time}, n::Integer, ϵ::R, tmpd::Vector{R}, primitive1!, primitive2!) where {R}
    return LanczosTimeTmp(Int(n), R(log(n)), Matrix{R}(undef, n, n), Vector{R}(undef, n -1), Vector{R}(undef, n -1),
        Vector{R}(undef, n), tmpd, primitive1!, primitive2!)
end

function approx_min_evec(t, idx, tmp::LanczosTimeTmp{R}) where {R}
    vmat, ρ, ω, v, tmpd, primitive1!, primitive2! = tmp.vmat, tmp.ρ, tmp.ω, tmp.v, tmp.tmpd, tmp.primitive1!, tmp.primitive2!
    @views @inbounds begin
        # v₁ ← randn(n, 1)
        # v₁ ← v₁ / ‖ v₁ ‖
        normalize!(Random.randn!(vmat[:, 1]))
        local i
        # for i ← 1, 2, 3, ..., min{q, n -1} do
        for outer i in 1:min(ceil(Int, t^R(1//4) * tmp.logn), tmp.n -1)
            vᵢ = vmat[:, i]
            vᵢ₊₁ = vmat[:, i+1]
            # ωᵢ ← Re(vᵢ* (M vᵢ)). So we first need M vᵢ, which is C + A*(y + β(z - b))
            primitive1!(vᵢ₊₁, vᵢ, idx, true, false)
            primitive2!(vᵢ₊₁, vᵢ, tmpd, idx, true, true)
            ω[i] = real(dot(vᵢ, vᵢ₊₁))
            # vᵢ₊₁ ← M vᵢ - ωᵢ vᵢ - ρᵢ₋₁ vᵢ₋₁ where ρ₀ v₀ = 0
            if i == 1
                vᵢ₊₁ .-= ω[i] .* vᵢ
            else
                vᵢ₊₁ .-= ω[i] .* vᵢ .+ ρ[i-1] .* vmat[:, i-1]
            end
            # ρᵢ ← ‖ vᵢ ‖ (should be vᵢ₊₁)
            ρ[i] = norm(vᵢ₊₁)
            # if ρᵢ = 0 then break
            ρ[i] < 1e-8 && break
            vᵢ₊₁ ./= ρ[i]
        end
        # T ← tridiag(ρ[1:i-1], ω[1:i], ρ[1:i-1])
        # [ξ, u] ← MinEvec(T)
        w, Z = LAPACK.stegr!('V', 'I', ω[1:i], ρ[1:i-1], 0.0, 0.0, 1, 1) # TODO: call LAPACK ourselves, so we avoid allocations
        ξ = w[1]
        # v ← ∑ⱼⁱ uⱼ vⱼ
        mul!(v, vmat[:, 1:i], Z[:, 1])
    end
    nrm = norm(v)
    return ξ * nrm, lmul!(inv(nrm), v)
end

struct LanczosSpaceTmp{R,P1,P2}
    n::Int
    logn::R
    v₁::Vector{R}
    v₂::Vector{R}
    v₃::Vector{R}
    ρ::Vector{R}
    ω::Vector{R}
    v::Vector{R}
    tmpd::Vector{R}
    primitive1!::P1
    primitive2!::P2
end

function setup_approx_min_evec(::Val{:lanczos_space}, n::Integer, ϵ::R, tmpd::Vector{R}, primitive1!, primitive2!) where {R}
    return LanczosSpaceTmp(Int(n), R(log(n)), Vector{R}(undef, n), Vector{R}(undef, n), Vector{R}(undef, n),
        Vector{R}(undef, n -1), Vector{R}(undef, n -1), Vector{R}(undef, n), tmpd, primitive1!, primitive2!)
end

function approx_min_evec(t, idx, tmp::LanczosSpaceTmp{R}) where {R}
    vᵢ₋₁, vᵢ, vᵢ₊₁, ρ, ω, v, tmpd, primitive1!, primitive2! = tmp.v₁, tmp.v₂, tmp.v₃, tmp.ρ, tmp.ω, tmp.v, tmp.tmpd,
        tmp.primitive1!, tmp.primitive2!
    @views @inbounds begin
        # v₁ ← randn(n, 1)
        # v₁ ← v₁ / ‖ v₁ ‖
        normalize!(Random.randn!(v))
        copyto!(vᵢ, v)
        local i
        # for i ← 1, 2, 3, ..., min{q, n -1} do
        for outer i in 1:min(ceil(Int, t^R(1//4) * tmp.logn), tmp.n -1)
            # ωᵢ ← Re(vᵢ* (M vᵢ)). So we first need M vᵢ, which is C + A*(y + β(z - b))
            primitive1!(vᵢ₊₁, vᵢ, idx, true, false)
            primitive2!(vᵢ₊₁, vᵢ, tmpd, idx, true, true)
            ω[i] = real(dot(vᵢ, vᵢ₊₁))
            # vᵢ₊₁ ← M vᵢ - ωᵢ vᵢ - ρᵢ₋₁ vᵢ₋₁ where ρ₀ v₀ = 0
            if i == 1
                vᵢ₊₁ .-= ω[i] .* vᵢ
            else
                vᵢ₊₁ .-= ω[i] .* vᵢ .+ ρ[i-1] .* vᵢ₋₁
            end
            # ρᵢ ← ‖ vᵢ ‖ (should be vᵢ₊₁)
            ρ[i] = norm(vᵢ₊₁)
            # if ρᵢ = 0 then break
            ρ[i] < 1e-8 && break
            vᵢ₊₁ ./= ρ[i]
            vᵢ₋₁, vᵢ, vᵢ₊₁ = vᵢ, vᵢ₊₁, vᵢ₋₁ # just reference cross-swapping
        end
        # T ← tridiag(ρ[1:i-1], ω[1:i], ρ[1:i-1])
        # [ξ, u] ← MinEvec(T)
        w, Z = LAPACK.stegr!('V', 'I', ω[1:i], ρ[1:i-1], 0.0, 0.0, 1, 1) # TODO: call LAPACK ourselves, so we avoid allocations
        ξ = w[1]
        # as we don't have the vectors available, we must redo all of the computation. We stored our initial random vector in
        # v, which we now flip back to vᵢ. Note that LAPACK overwrites both ω and ρ, so we also have to re-compute those.
        vᵢ, v = v, vᵢ
        local ωᵢ, ρᵢ
        for i in 1:i
            primitive1!(vᵢ₊₁, vᵢ, idx, true, false)
            primitive2!(vᵢ₊₁, vᵢ, tmpd, idx, true, true)
            ωᵢ = real(dot(vᵢ, vᵢ₊₁)) # no need to do some index calculations, just a simple variable
            if i == 1
                vᵢ₊₁ .-= ωᵢ .* vᵢ
                v .= Z[1, 1] * vᵢ
            else
                vᵢ₊₁ .-= ωᵢ .* vᵢ .+ ρᵢ .* vᵢ₋₁ # indeed, this is ρᵢ₋₁ now, but we don't introduce an alias
                v .+= Z[i, 1] * vᵢ
            end
            ρᵢ = norm(vᵢ₊₁)
            ρᵢ < 1e-8 && break
            vᵢ₊₁ ./= ρᵢ
            vᵢ₋₁, vᵢ, vᵢ₊₁ = vᵢ, vᵢ₊₁, vᵢ₋₁ # just reference cross-swapping
        end
    end
    nrm = norm(v)
    return ξ * nrm, lmul!(inv(nrm), v)
end

struct MatrixFreeOperator{T,P1}
    size::Tuple{Int,Int}
    primitive!::P1
end

MatrixFreeOperator{R}(size::Tuple{Int,Int}, primitive!::P1) where {R,P1} = MatrixFreeOperator{R,P1}(size, primitive!)

Base.@propagate_inbounds LinearAlgebra.mul!(C, AorB::MatrixFreeOperator, X, α, β) = AorB.primitive!(C, X, α, β)
Base.size(A::MatrixFreeOperator, i::Int) = A.size[i]
Base.eltype(::MatrixFreeOperator{T}) where {T} = T

struct LOBPCGFastTmp{R,I}
    n::Int
    logn::R
    iterator::I
    idx::Base.RefValue{Int}
    v::Vector{R}
    tmpd::Vector{R}
end

function setup_approx_min_evec(::Val{:lobpcg_fast}, n::Integer, ϵ::R, tmpd::Vector{R}, primitive1!, primitive2!) where {R}
    v = Vector{R}(undef, n)
    idx = Ref{Int}()
    return LOBPCGFastTmp(
        Int(n),
        R(log(n)),
        IterativeSolvers.LOBPCGIterator(
            MatrixFreeOperator{R}((n, n), (out, in, alpha, beta) -> begin
                primitive1!(out, in, idx[], alpha, beta)
                primitive2!(out, in, tmpd, idx[], alpha, true)
                nothing
            end), false, reshape(v, :, 1)
        ),
        idx,
        v,
        tmpd
    )
end

function approx_min_evec(t, idx, tmp::LOBPCGFastTmp{R}) where {R}
    normalize!(Random.randn!(tmp.v))
    tmp.idx[] = idx
    ξ = IterativeSolvers.lobpcg!(tmp.iterator, tol=1e-8, maxiter=min(ceil(Int, t^R(1//4) * tmp.logn), tmp.n -1)).λ[1]
    return ξ, tmp.v
end

struct LOBPCGAccurateTmp{R,I}
    ϵ::R
    iterator::I
    idx::Base.RefValue{Int}
    v::Vector{R}
    tmpd::Vector{R}
end

function setup_approx_min_evec(::Val{:lobpcg_accurate}, n::Integer, ϵ, tmpd::Vector{R}, primitive1!, primitive2!) where {R}
    v = Vector{R}(undef, n)
    idx = Ref{Int}()
    return LOBPCGAccurateTmp(
        iszero(ϵ) ? IterativeSolvers.default_tolerance(R) : ϵ/100,
        IterativeSolvers.LOBPCGIterator(
            MatrixFreeOperator{R}((n, n), (out, in, alpha, beta) -> begin
                primitive1!(out, in, idx[], alpha, beta)
                primitive2!(out, in, tmpd, idx[], alpha, true)
                nothing
            end), false, reshape(v, :, 1)
        ),
        idx,
        v,
        tmpd
    )
end

function approx_min_evec(t, idx, tmp::LOBPCGAccurateTmp)
    normalize!(Random.randn!(tmp.v))
    tmp.idx[] = idx
    ξ = IterativeSolvers.lobpcg!(tmp.iterator, tol=tmp.ϵ).λ[1]
    return ξ, tmp.v
end

struct NystromSketch{T}
    Ω::Vector{Matrix{T}}
    S::Vector{Matrix{T}}
    tmp::Vector{Vector{T}}

    function NystromSketch{T}(n::AbstractVector{<:Integer}, R::AbstractVector{<:Integer}) where {T}
        @assert(length(n) == length(R))
        tmp = Vector{T}(undef, maximum(R))
        @inbounds return new(
            [randn(T, nᵢ, Rᵢ) for (nᵢ, Rᵢ) in zip(n, R)],
            [zeros(T, nᵢ, Rᵢ) for (nᵢ, Rᵢ) in zip(n, R)],
            [@view(tmp[1:Rᵢ]) for Rᵢ in R]
        )
    end
end

# Note: This is a sequential implementation (we also share temporaries). The bottleneck is the memory access, so we don't try
# any parallelization, although every iteration is (provided the temporaries are separated) independent from the other.
rank_one_update!(sketch::NystromSketch{T}, v::AbstractVector{T}, i::Integer, η::T) where {T} =
    # S ← (1 - η) S + η v (v* Ω)
    mul!(sketch.S[i], v, transpose(mul!(sketch.tmp[i], transpose(sketch.Ω[i]), conj(v))), η, 1 - η)
    # TODO: benchmark whether scaling + rank-1 update is better than level 3 function.

maxcolnorm(M::AbstractMatrix{T}) where {T} = maximum(norm, _eachcol(M); init=zero(T))::T

function reconstruct(sketch::NystromSketch{T}) where {T}
    n = size.(sketch.Ω, 1)
    # σ ← √n eps(norm(S))
    σ = (sqrt.(n) .* eps(T) .* maxcolnorm.(sketch.S))::Vector{T}
    # Sσ ← S + σ Ω
    Sₛ = (sketch.S .+ σ .* sketch.Ω)::Vector{Matrix{T}}
    # L ← chol(Ω* Sσ).
    L = cholesky!.(Positive, adjoint.(sketch.Ω) .* Sₛ)::Vector{Cholesky{T,Matrix{T}}}
    # [U, Σ, ~] ≤ svd(Sσ / L)
    svds = svd!.(rdiv!.(Sₛ, getproperty.(L, :U)))::Vector{SVD{T,T,Matrix{T},Vector{T}}}
    # Λ ← max{0, Σ² - σ I}
    Λ = [max.(zero(T), svdᵢ.S .^ 2 .- σᵢ)::Vector{T} for (svdᵢ, σᵢ) in zip(svds, σ)]

    return getproperty.(svds, :U)::Vector{Matrix{T}}, Λ
end

end
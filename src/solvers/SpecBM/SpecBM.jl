export SpecBMResult, specbm_primal

# This is an implementation of the SpecPM primal solver, https://arxiv.org/abs/2307.07651v1 with a reference implementation on
# https://github.com/soc-ucsd/SpecBM.git, tightly integrated with the PolynomialOptimization framework
const VecView{R,cont} = SubArray{R,1,Vector{R},Tuple{UnitRange{Int}},cont}
const MatView{R,cont} = SubArray{R,2,Matrix{R},Tuple{UnitRange{Int},UnitRange{Int}},cont}
const DiagView{R} = SubArray{R,1,Base.ReshapedArray{R,1,MatView{R,false},Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int}}},Tuple{StepRange{Int,Int}},false}

struct SpecBMData{R,PType,AType,AtType,AVType,APVType,BType,CType,CVType}
    psds::PType
    r::Vector{Int}
    ϵ::R

    Ω::Vector{R}
    w_psd::Vector{R}
    P_psds::Vector{Matrix{R}}

    A::AType
    At::AtType
    a_free::AVType
    a_psd::AVType
    a_psds::APVType
    b::BType
    c::CType
    c_free::CVType
    c_psd::CVType
    C_psds::Vector{PackedMatrix{R,CVType,:LS}}
    ω_free::VecView{R,true}
    ω_psd::VecView{R,true}
    Ω_psds::Vector{PackedMatrix{R,VecView{R,true},:LS}}
    W_psds::Vector{PackedMatrix{R,VecView{R,true},:LS}}

    function SpecBMData(num_vars::Integer, num_frees::Integer, psds::AbstractVector{<:Integer}, r::Vector{Int}, ϵ::R,
        A::AbstractMatrix{R}, At::AbstractMatrix{R}, b::AbstractVector{R}, c::AbstractVector{R}) where {R}
        @inbounds @views begin
            @assert(length(psds) == length(r))
            num_psdvars = sum(packedsize, psds, init=0)
            @assert(num_frees + num_psdvars == num_vars)
            num_psds = length(r)
            # allocated problem data
            Ω = zeros(R, num_vars)
            w_psd = zeros(R, num_psdvars)
            W_psds = Vector{PackedMatrix{R,typeof(w_psd[begin:end]),:LS}}(undef, num_psds)
            P_psds = Vector{Matrix{R}}(undef, num_psds)
            # views into existing data
            a_free = A[:, 1:num_frees]
            a_psd = A[:, num_frees+1:end]
            a_psds = Vector{typeof(A[:, begin:end])}(undef, num_psds)
            c_free = c[1:num_frees]
            c_psd = c[num_frees+1:end]
            C_psds = Vector{PackedMatrix{R,typeof(c[begin:end]),:LS}}(undef, num_psds)
            ω_free = Ω[1:num_frees]
            ω_psd = Ω[num_frees+1:end]
            Ω_psds = Vector{PackedMatrix{R,typeof(Ω[begin:end]),:LS}}(undef, num_psds)
            i = num_frees +1
            for (j, (nⱼ, rⱼ)) in enumerate(zip(psds, r))
                # initialize all the data and connect the views appropriately
                dimⱼ = packedsize(nⱼ)
                Ω_psds[j] = Ωⱼ = PackedMatrix(nⱼ, Ω[i:i+dimⱼ-1], :LS)
                # An initial point Ω₀ ∈ 𝕊ⁿ.  As in the reference implementation, we take zero for the free variables and the
                # vectorized identity for the PSD variables.
                for k in PackedDiagonalIterator(Ωⱼ)
                    Ωⱼ[k] = one(R)
                end
                # Initialize W̄₀ ∈ 𝕊₊ⁿ with tr(W̄₀) = 1. As in the reference implementation, we take the (1,1) elementary matrix.
                # Note that the reference implementation only allows for a single block; we map this to multiple semidefinite
                # constraints not merely by mimicking a block-diagonal matrix, but taking the constraints into account
                # individually!
                W_psds[j] = Wⱼ = PackedMatrix(nⱼ, w_psd[i-num_frees:i-num_frees+dimⱼ-1], :LS)
                Wⱼ[1, 1] = one(R)
                # Compute P₀ ∈ ℝⁿˣʳ with columns being the top r orthonormal eigenvectors of -Ω₀. As Ω₀ is the identity, we can
                # do this explicitly.
                P_psds[j] = Pⱼ = zeros(R, nⱼ, rⱼ)
                for k in 1:rⱼ
                    Pⱼ[k, k] = one(R)
                end
                a_psds[j] = A[:, i:i+dimⱼ-1]
                C_psds[j] = PackedMatrix(nⱼ, c[i:i+dimⱼ-1], :LS)

                i += dimⱼ
            end
        end

        return new{R,typeof(psds),typeof(A),typeof(At),typeof(a_free),typeof(a_psds),typeof(b),typeof(c),typeof(c_free)}(
            psds, r, ϵ,
            Ω, w_psd, P_psds,
            A, At, a_free, a_psd, a_psds, b, c, c_free, c_psd, C_psds, ω_free, ω_psd, Ω_psds, W_psds
        )
    end
end

function Base.getproperty(d::SpecBMData, name::Symbol)
    name === :num_vars && return length(getfield(d, :Ω))
    name === :num_conds && return size(getfield(d, :A), 1)
    name === :num_frees && return size(getfield(d, :a_free), 2)
    name === :num_psds && return length(getfield(d, :psds))
    return getfield(d, name)
end
Base.propertynames(::SpecBMData) = (:num_vars, :num_conds, :num_frees, :num_psds, fieldnames(SpecBMData)...)

struct SpecBMMastersolverData{R}
    Xstar::Vector{R}
    sstar_psd::Vector{R}
    γstars::Vector{R}
    ystar::Vector{R}
    wstar_psd::Vector{R}

    xstar_free::VecView{R,true}
    xstar_psd::VecView{R,true}
    Xstar_psds::Vector{PackedMatrix{R,VecView{R,true},:LS}}
    Sstar_psds::Vector{PackedMatrix{R,VecView{R,true},:LS}}
    Wstar_psds::Vector{PackedMatrix{R,VecView{R,true},:LS}}

    function SpecBMMastersolverData(data::SpecBMData{R}) where {R}
        @inbounds begin
            num_psds = data.num_psds
            num_conds = data.num_conds
            num_frees = data.num_frees
            # allocated mastersolver output data
            Xstar = similar(data.Ω)
            sstar_psd = Vector{R}(undef, sum(packedsize, data.r, init=0))
            γstars = Vector{R}(undef, num_psds)
            ystar = Vector{R}(undef, num_conds)
            wstar_psd = similar(data.w_psd)
            # views into existing data
            xstar_free = @view(Xstar[1:num_frees])
            xstar_psd = @view(Xstar[num_frees+1:end])
            Xstar_psds = Vector{PackedMatrix{R,typeof(@view(Xstar[begin:end])),:LS}}(undef, num_psds)
            Sstar_psds = Vector{PackedMatrix{R,typeof(@view(sstar_psd[begin:end])),:LS}}(undef, num_psds)
            Wstar_psds = Vector{PackedMatrix{R,typeof(@view(wstar_psd[begin:end])),:LS}}(undef, num_psds)

            i_n = num_frees +1
            i_r = 1
            for (j, (nⱼ, rⱼ)) in enumerate(zip(data.psds, data.r))
                dimⱼ = packedsize(nⱼ)
                Xstar_psds[j] = PackedMatrix(nⱼ, @view(Xstar[i_n:i_n+dimⱼ-1]), :LS)
                Wstar_psds[j] = PackedMatrix(nⱼ, @view(wstar_psd[i_n-num_frees:i_n-num_frees+dimⱼ-1]), :LS)
                i_n += dimⱼ
                rdimⱼ = packedsize(rⱼ)
                Sstar_psds[j] = PackedMatrix(rⱼ, @view(sstar_psd[i_r:i_r+rdimⱼ-1]), :LS)
                i_r += rdimⱼ
            end
        end

        return new{R}(
            Xstar, sstar_psd, γstars, ystar, wstar_psd,
            xstar_free, xstar_psd, Xstar_psds, Sstar_psds, Wstar_psds
        )
    end
end

struct SpecBMCache{R,F,ACV,SS}
    # data for the actual minimization
    m₁::Vector{R}
    m₂::Vector{R}
    M::Symmetric{R,Matrix{R}}
    # views into the data
    M₁₁::MatView{R,false}
    M₂₁::MatView{R,false}
    M₂₂::MatView{R,false}
    # data/views for the preprocessing stage
    Pkrons::Vector{Matrix{R}}
    m₂s::Vector{VecView{R,true}}
    q₃::Vector{R}
    Q₁₁::DiagView{R} # diagonal of M₁₁
    Q₂₁s::Vector{SubArray{R,1,Matrix{R},Tuple{UnitRange{Int},Int},true}} # block-diagonal
    Q₂₂::DiagView{R} # diagonal of M₂₂
    Q₃₁::Matrix{R}
    Q₃₂::Matrix{R}
    Q₃₂s::Vector{SubArray{R,2,Matrix{R},Tuple{Base.Slice{Base.OneTo{Int}},UnitRange{Int}},true}}
    Q₃₃inv::F
    # some precomputed data
    Σr::Int
    twoAc::ACV
    # caches for eigendecomposition
    eigens::Vector{Tuple{Eigen{R,R,Matrix{R},Vector{R}},Vector{R},Vector{BLAS.BlasInt},Vector{BLAS.BlasInt},Matrix{R}}}
    # and one temporary in various forms (shared memory!)
    tmp::Vector{R}
    # finally the subsolver
    subsolver::SS

    function SpecBMCache(data::SpecBMData{R}, AAt, subsolver, ρ, r_current) where {R}
        @inbounds begin
            rdims = packedsize.(data.r)
            Σr = sum(rdims, init=0)
            num_psds = data.num_psds
            num_conds = data.num_conds
            # allocated minimization data
            m₁ = Vector{R}(undef, num_psds)
            m₂ = Vector{R}(undef, Σr)
            M = Matrix{R}(undef, num_psds + Σr, num_psds + Σr)
            # views into the data
            M₁₁ = @view(M[1:num_psds, 1:num_psds])
            M₂₁ = @view(M[num_psds+1:end, 1:num_psds])
            M₂₂ = @view(M[num_psds+1:end, num_psds+1:end])

            # data/views for the preprocessing stage
            Pkrons = Vector{Matrix{R}}(undef, num_psds)
            m₂s = Vector{typeof(@view(m₂[begin:end]))}(undef, num_psds)
            q₃ = Vector{R}(undef, num_conds)
            Q₁₁ = @view(M₁₁[begin:num_psds+1:end])
            Q₂₁s = Vector{typeof(@view(M₂₁[begin:end, begin]))}(undef, num_psds)
            Q₂₂ = @view(M₂₂[begin:Σr+1:end])
            Q₃₁ = Matrix{R}(undef, num_conds, num_psds)
            Q₃₂ = Matrix{R}(undef, num_conds, Σr)
            Q₃₂s = Vector{typeof(@view(Q₃₂[:, begin:end]))}(undef, num_psds)
            Q₃₃inv = try EfficientCholmod(ldlt(AAt)) catch; qr(AAt) end
            twoAc = rmul!(data.A * data.c, R(2)) # typically, A and c are sparse, so the * implementation is the best
            eigens = Vector{Tuple{Eigen{R,R,Matrix{R},Vector{R}},Vector{R},Vector{BLAS.BlasInt},Vector{BLAS.BlasInt},Matrix{R}}}(undef, num_psds)
            tmp = Vector{R}(undef, max(num_conds * max(num_psds, Σr), maximum(data.r, init=0)^2, maximum(data.psds, init=0)^2))
            i = 1
            for (j, (nⱼ, rⱼ, rdimⱼ, r_currentⱼ)) in enumerate(zip(data.psds, data.r, rdims, r_current))
                Pkrons[j] = Matrix{R}(undef, packedsize(nⱼ), rdimⱼ)
                m₂s[j] = @view(m₂[i:i+rdimⱼ-1])
                Q₂₁s[j] = @view(M₂₁[i:i+rdimⱼ-1, j])
                Q₃₂s[j] = @view(Q₃₂[:, i:i+rdimⱼ-1])
                i += rdimⱼ
                eigens[j] = ( # we need nⱼ buffer space for the eigenvalues
                    Eigen(Vector{R}(undef, nⱼ), Matrix{R}(undef, nⱼ, min(r_currentⱼ, nⱼ))),
                    Vector{R}(undef, max(8nⱼ, 1 + 6rⱼ + rⱼ^2)),
                    Vector{BLAS.BlasInt}(undef,  max(5nⱼ, 3 + 5rⱼ)),
                    Vector{BLAS.BlasInt}(undef, nⱼ),
                    Matrix{R}(undef, rⱼ, rⱼ)
                )
                # this is not excessive - LAPACK requires nⱼ buffer space for the eigenvalues even if less are requested
                # while if r_currentⱼ == nⱼ, we will call spevd! instead of spevx! which has different workspace rules. But we
                # also need to find the minimum eigenvalue of Ωⱼ, for which we always call spevx!, and we always need the full
                # eigendecomposition of Sⱼ with spevd!.
            end
            ss = specbm_setup_primal_subsolver(Val(subsolver), num_psds, data.r, rdims, Σr, ρ)
        end

        return new{R,typeof(Q₃₃inv),typeof(twoAc),typeof(ss)}(
            m₁, m₂, Symmetric(M, :L),
            M₁₁, M₂₁, M₂₂,
            Pkrons, m₂s, q₃, Q₁₁, Q₂₁s, Q₂₂, Q₃₁, Q₃₂, Q₃₂s, Q₃₃inv,
            Σr, twoAc,
            eigens,
            tmp,
            ss
        )
    end
end

gettmp(c::SpecBMCache, sizes...) = reshape(@view(c.tmp[1:*(sizes...)]), sizes...)

function Base.getproperty(c::SpecBMCache, name::Symbol)
    name === :q₁ && return getfield(c, :m₁)
    name === :q₂s && return getfield(c, :m₂s)
    return getfield(c, name)
end
Base.propertynames(::SpecBMCache) = (:q₁, :q₂s, fieldnames(SpecBMCache)...)

"""
    SpecBMResult

Contains the result of a SpecBM run

# Fields
- `status::Symbol`: one of `:Optimal`, `:IterationLimit`, `:SlowProgress`
- `objective::R`: the objective value
- `x::Vector{R}`: the optimal vector of primal variables: first, `num_frees` free variables, then all scaled vectorized lower
  triangles of the PSD variables
- `y::Vector{R}`: the optimal vector of dual variables, one for each constraint
- `iterations::Int`: the number of iterations until the given status was reached
- `quality::R`: the optimality quantifier that is compared against `ϵ` to determine convergence, which is determined by the
  maximum of the relative quantities below and the negative primal infeasibility.
- `primal_infeas::R`
- `dual_infeas::R`
- `gap::R`
- `rel_accuracy::R`
- `rel_primal_infeas::R`
- `rel_dual_infeas::R`
- `rel_gap`
"""
struct SpecBMResult{R}
    status::Symbol
    objective::R
    x::Vector{R}
    y::Vector{R}
    iterations::Int
    quality::R
    primal_infeas::R
    dual_infeas::R
    gap::R
    rel_accuracy::R
    rel_primal_infeas::R
    rel_dual_infeas::R
    rel_gap::R
end

@doc raw"""
    specbm_primal(A, b, c; num_frees=missing, psds, ρ, r_past, r_current, ϵ=1e-4, β=0.1, maxiter=10000, maxnodescent=15,
        adaptiveρ=false, α=1., adaptiveα=true, αmin=1e-5, αmax=1000., ml=0.001, mu=min(1.5β, 1), Nmin=10, verbose=false,
        step=20, offset=0, At=transpose(A), AAt=A*At, subsolver=:Mosek, callback=(data, mastersolver_data)->nothing)

Solves the minimization problem
```math
    \min_x \{ ⟨c, x⟩ : A x = b, x = (x_{\mathrm{free}}, \operatorname{svec}(X_1), \dotsc), X_i ⪰ 0,
              \sum_i \operatorname{tr}(X_i) ≤ ρ \} + \mathit{offset}
```
where the vector ``x`` contains `num_frees` free variables, followed by the vectorized and scaled lower triangles of PSD
matrices ``X_i`` that have side dimensions given in `psds`. _Scaled_ here means that the off-diagonal elements must be
multiplied by ``\sqrt2`` when going from the matrix to its vectorization, so that scalar products are preserved. This
corresponds to the `:LS` format of a [`PackedMatrix`](@ref).

# Arguments
## Problem formulation
- `A::AbstractMatrix{R}`: a sparse or dense matrix
- `At::AbstractMatrix{R}`: the transpose of `A`. If omitted, `transpose(A)` is used instead. However, if the transpose is
  already known in explicit form (in particular, as another `SparseMatrixCSC`), some operations can be carried out faster.
- `AAt::AbstractMatrix{R}`: the product `A*At`, which is also calculated automatically, but can be given if it is already
  known.
- `b::AbstractVector{R}`: a dense or sparse vector
- `c::AbstractVector{R}`: a dense or sparse vector
- `offset::Real`: an offset that is added to the objective
- `num_frees`: the number of free variables in the problem. The first `num_frees` entries in ``x`` will be free. If this value
  is omitted, it is automatically calculated based on the dimensions of `A` and `psds`.
- `psds::AbstractVector{<:Integer}`: a vector that, for each semidefinite matrix in the problem, specifies its side dimension.
  A side dimension of ``n`` will affect ``\frac{n(n +1)}{2}`` variables.
- `ρ::Real`: an upper bound on the total trace in the problem. Note that by setting `adaptiveρ=true`, this bound will
  effectively be removed by dynamically growing as necessary. In this case, the value specified here is the initial value.
- `adaptiveρ::Bool`: effectively sets ``\rho \to \infty``; note that an initial `ρ` still has to be provided.
## Spectral bundle parameters
- `r_past::Integer`: the number of past eigenvectors to keep, must be nonnegative
- `r_current::Integer`: the number of current eigenvectors to keep, must be positive
- `β::Real`: A step is recognized as a descent step if the decrease in the objective value is at least a factor
  ``\beta \in (0, 1)`` smaller than the decrease predicted by the model.
- `α::Real`: the regularization parameter for the augmented Lagrangian; must be positive
- `adaptiveα::Bool=true`: enables adaptive updating of `α` depending on the following five parameters, as described in
  [Liao et al](https://doi.org/10.48550/arXiv.2307.07651).
- `αmin::Real`: lower bound for the adaptive algorithm that `α` may not exceed
- `αmax::Real`: upper bound for the adaptive algorithm that `α` may not exceed
- `ml::Real`: `α` is doubled if the decrease in the objective value is at least a factor ``m_{\mathrm l} \in (0, \beta)``
  larger than predicted by the model, provided no descent step was recognized for at least `Nmin` iterations.
- `mu::Real`: `α` is halved if the decrease in the objective value is at least a factor ``m_{\mathrm u} > \beta`` smaller than
  predicted by the model.
- `Nmin::Integer`: minimum number of no-descent-steps before `ml` becomes relevant
- `subsolver::Symbol`: subsolver to solve the quadratic semidefinite subproblem in every iteration of SpecBM. Currently,
  `:Hypatia` and `:Mosek` are supported; however, note that Mosek will require at least version 10.1.11 (better 10.1.13 to
  avoid some rare crashes).
## Termination criteria
- `ϵ::Real`: minimum quality of the result in order for the algorithm to terminate successfully (status `:Optimal`)
- `maxiter::Integer`: maximum number of iterations before the algorithm terminates anyway (status `:IterationLimit`). Must be
  at least `2`.
- `maxnodescent::Integer`: maximum number of consecutive iterations that may report no descent step before the algorithm
  terminates (status `:SlowProgress`). Must be positive or zero to disable this check.
## Logging
- `verbose::Bool`: print the status every `step` iterations. Note that the first (incomplete) iteration will never be printed.
- `step::Integer`: skip a number of iterations and only print every `step`th.
## Advanced solver interaction
- `callback::Function`: a callback that is called with the last problem data (type `SpecBMData`) and the last mastersolver data
  (type `SpecBMMastersolverData`) before the mastersolver is called anew. Changes to the structures may be made.

See also [`SpecBMResult`](@ref).
"""
function specbm_primal(A::AbstractMatrix{R}, b::AbstractVector{R}, c::AbstractVector{R};
    num_frees::Union{Missing,Integer}=missing, psds::AbstractVector{<:Integer},
    ρ::Real, r_past::Union{<:AbstractVector{<:Integer},<:Integer}, r_current::Union{<:AbstractVector{<:Integer},<:Integer},
    ϵ::Real=R(1e-4), β::Real=R(0.1), maxiter::Integer=10000, maxnodescent::Integer=15, adaptiveρ::Bool=false,
    α::Real=R(1.), adaptiveα::Bool=true, αmin::Real=R(1e-5), αmax::Real=R(1000.),
    ml::Real=R(0.001), mr::Real=min(R(1.5) * β, 1), Nmin::Integer=10,
    verbose::Bool=true, step::Integer=20, offset::R=zero(R),
    At::Union{Missing,AbstractMatrix{R}}=missing, AAt::Union{Missing,AbstractMatrix{R}}=missing,
    subsolver::Symbol=:Mosek, callback::Function=(data, mastersolver_data) -> nothing) where {R<:AbstractFloat}
    #region Input validation
    subsolver ∈ (:Mosek, :Hypatia) || error("Unsupported subsolver ", subsolver)
    # Problem data A₁, ..., Aₘ, C ∈ 𝕊ⁿ, b ∈ ℝⁿ. Here, we also allow for free variables, as in the reference implementation.
    # We do not store the matrices A directly, but instead interpret all PSD variables by their scaled vectorized upper
    # triangle (contrary to the reference implementation, which uses vectorized full storage). Therefore, A contains the
    # (row-wise) stacked vectorized matrices with off-diagonals scaled by √2 and C is also a vector similarly scaled. All free
    # variables come before the PSD variables.
    num_conds, num_vars = size(A)
    (num_conds == length(b) && num_vars == length(c)) || error("Incompatible dimensions")
    all(j -> j > 0, psds) || error("PSD dimensions must be positive")
    if ismissing(num_frees)
        num_frees = num_vars - sum(packedsize, psds, init=0)
        num_frees ≥ 0 || error("Incompatible dimensions")
    elseif num_frees < 0
        error("Number of free variables must be nonnegative")
    elseif sum(packedsize, psds, init=0) + num_frees != num_vars
        error("Incompatible dimensions")
    end
    num_psds = length(psds)
    if isa(r_current, Integer)
        r_current ≥ 0 || error("r_current must be positive")
        r_current = min.(r_current, psds)
    elseif length(r_current) != num_psds
        error("Number of r_current must be the same as number of psd constraints")
    else
        all(x -> x ≥ 1, r_current) || error("r_current must be positive")
        all(splat(≤), zip(r_current, psds)) || error("No r_current must not exceed its associated dimension")
    end
    if isa(r_past, Integer)
        r_past ≥ 0 || error("r_past must be nonnegative")
        r_past = min.(fill(r_past, num_psds), psds .- r_current) # which is guaranteed to be nonnegative
    elseif length(r_past) != num_psds
        error("Number of r_past must be the same as number of psd constraints")
    else
        all(x -> x ≥ 0, r_past) || error("r_past must be nonnegative")
        all((r_currentⱼ, r_pastⱼ, dimⱼ) -> r_currentⱼ + r_pastⱼ ≤ dimⱼ) ||
            error("r_past + r_current must not exceed the associated dimension")
    end
    # Parameters rₚ ≥ 0, r_c ≥ 1, α > 0, β ∈ (0, 1), ϵ ≥ 0, tₘₐₓ ≥ 1
    α > 0 || error("α must be positive")
    0 < β < 1 || error("β must be in (0, 1)")
    ϵ ≥ 0 || error("ϵ must be nonnegative")
    maxiter > 1 || error("maxiter must be larger than 1")
    maxnodescent ≥ 0 || error("maxnodescent must be nonnegative")
    # Adaptive parameters mᵣ > β, 0 < mₗ < β
    if adaptiveα
        mr > β || error("mr must be larger than β")
        0 < ml < β || error("ml must be in (0, β)")
        0 < Nmin || error("Nmin must be positive")
        iszero(maxnodescent) || maxnodescent ≥ Nmin || error("maxnodescend must not be smaller than Nmin")
        α = inv(R(2))
    end
    if ismissing(At)
        At = transpose(A) # it would be best if A already was a transpose(At), as we need slices of rows in A
    end
    if ismissing(AAt)
        AAt = A * At
    end
    step ≥ 1 || error("step must be positive")
    #endregion

    @verbose_info("SpecBM Primal Solver with parameters ρ = $ρ, r_past = $r_past, r_current = $r_current, ϵ = $ϵ, β = $β, $α ",
        adaptiveα ? "∈ [$αmin, $αmax], ml = $ml, mr = $mr" : "= $α", ", subsolver = $subsolver")
    @verbose_info("Iteration | Primal objective | Primal infeas | Dual infeas | Duality gap | Rel. accuracy | Rel. primal inf. | Rel. dual inf. |    Rel. gap | Descent step | Consecutive null steps",
        adaptiveρ ? " | Dual trace" : "")

    invnormbplus1 = inv(norm(b) + one(R))
    invnormcplus1 = inv(norm(c) + one(R))

    data = SpecBMData(num_vars, num_frees, psds, Int.(r_past .+ r_current), ϵ, A, At, b, c)
    mastersolver = SpecBMMastersolverData(data)
    cache = SpecBMCache(data, AAt, subsolver, ρ, r_current)

    # We need some additional variables for the adaptive strategy, following the naming in the reference implementation
    # (in the paper, the number of consecutive null steps N_c is used instead).
    null_count = 0
    has_descended = true
    changed_ρ = false
    status = :IterationLimit

    # 2: for t = 0, ..., tₘₐₓ do [we fix this to 1:maxiter]
    local FΩ, t, quality, primal_infeas, dual_infeas, gap, rel_accuracy, rel_primal_infeas, rel_dual_infeas, rel_gap
    for outer t in 1:maxiter
        t > 2 && callback(data, mastersolver)
        # 3: solve (24) to obtain Xₜ₊₁*, γₜ*, Sₜ*
        # combined with
        # 4: form the iterate Wₜ* in (28) and dual iterate yₜ* in (29)
        dfeasi, dfeasi_psd, dfeasi_free, gap = direction_qp_primal_free!(mastersolver, data, !isone(t), α, cache)
        # We also calculate some quality criteria here
        dual_infeas = max(dfeasi_free, dfeasi_psd)
        rel_dual_infeas = sqrt(dfeasi * invnormcplus1)
        if has_descended
            rel_primal_infeas = let tmp=gettmp(cache, length(b))
                copyto!(tmp, b)
                mul!(tmp, A, data.Ω, true, -one(R))
                norm(tmp) * invnormbplus1
            end
            # else we no not need to recompute this, the value from the last iteration is still valid
        end
        # 5: if t = 0 and A(Ωₜ) ≠ b then
        if isone(t) && rel_primal_infeas > ϵ # note: reference implementation does not check A(Ωₜ) ≠ b
            copyto!(data.Ω, mastersolver.Xstar)
            # we need the eigendecomposition for later in every case
            for (j, ((ev, work, iwork, ifail, _), Xstarⱼ)) in enumerate(zip(cache.eigens, mastersolver.Xstar_psds))
                if ==(size(ev.vectors)...)
                    eigen!(Xstarⱼ, ev.values, ev.vectors, work)
                else
                    @inbounds eigen!(Xstarⱼ, 1:r_current[j], ev.values, ev.vectors, work, iwork, ifail)
                end
            end
        # 7: else
        else
            # 8: if (25) holds then
            # (25): β( F(Ωₜ) - ̂F_{Wₜ, Pₜ}(Xₜ₊₁*)) ≤ F(Ωₜ) - F(Xₜ₊₁*)
            # where (20): F(X) := ⟨C, X⟩ - ρ min(λₘᵢₙ(X), 0)
            if has_descended || changed_ρ
                Σ = zero(R)
                for ((ev, work, iwork, ifail, _), Ωⱼ) in zip(cache.eigens, data.Ω_psds)
                    Ωcopy = PackedMatrix(LinearAlgebra.checksquare(Ωⱼ), gettmp(cache, length(Ωⱼ)),
                        PackedMatrices.packed_format(Ωⱼ))
                    copyto!(Ωcopy, Ωⱼ)
                    Σ += min(eigmin!(Ωcopy, ev.values, ev.vectors, work, iwork, ifail), zero(R))
                end
                FΩ = dot(data.c, data.Ω) - ρ * Σ
                # else we do not need to recalculate this, it did not change from the previous iteration
            end
            cXstar = dot(data.c, mastersolver.Xstar)
            Fmodel = cXstar - dot(mastersolver.wstar_psd, mastersolver.xstar_psd)
            Σ = zero(R)
            for (j, ((ev, work, iwork, ifail), Xstarⱼ)) in enumerate(zip(cache.eigens, mastersolver.Xstar_psds))
                Xcopy = PackedMatrix(LinearAlgebra.checksquare(Xstarⱼ), gettmp(cache, length(Xstarⱼ)),
                    PackedMatrices.packed_format(Xstarⱼ))
                copyto!(Xcopy, Xstarⱼ)
                if ==(size(ev.vectors)...)
                    eigen!(Xcopy, ev.values, ev.vectors, work, iwork)
                else
                    @inbounds eigen!(Xcopy, 1:r_current[j], ev.values, ev.vectors, work, iwork, ifail)
                end
                Σ += min(first(ev.values), zero(R))
            end
            FXstar = cXstar - ρ * Σ
            estimated_drop = FΩ - Fmodel
            cost_drop = FΩ - FXstar
            if (has_descended = (β * estimated_drop ≤ cost_drop))
                # 9: set primal iterate Ωₜ₊₁ = Xₜ₊₁*
                copyto!(data.Ω, mastersolver.Xstar)
                # 6.1.1. Adaptive strategy (can only be lower case due to mₗ < β < mᵣ)
                if adaptiveα
                    if mr * estimated_drop ≤ cost_drop
                        α = max(α / 2, αmin)
                    end
                    null_count = 0
                end
            # 10: else
            else
                # 11: set primal iterate Ωₜ₊₁ = Ωₜ (no-op)
                # 6.1.1. Adaptive strategy (can only be upper case)
                if adaptiveα
                    null_count += 1
                    if null_count ≥ Nmin && ml * estimated_drop ≥ cost_drop
                        α = min(2α, αmax)
                        null_count = 0
                    end
                end
            # 12: end if
            end
            rel_accuracy = estimated_drop / (abs(FΩ) + one(R))
        # 13: end if
        end
        rel_gap = gap / (one(R) + abs(dot(data.c, data.Ω)) + abs(dot(data.b, mastersolver.ystar))) # now Ω is corrected
        # 14: compute Pₜ₊₁ as (26), and Wₜ₊₁ as (27)
        # (26): Pₜ₊₁ = orth([Vₜ; Pₜ Q₁])
        # where Vₜ: top r_c ≥ 1 eigenvectors of -Xₜ₊₁*
        # and S* = [Q₁ Q₂] * Diagonal(Σ₁, Σ₂) * [Q₁; Q₂] with division in (rₚ, r - rₚ)
        # (27): Wₜ₊₁ = 1/(γ* + tr(Σ₂)) * (γ* Wₜ + Pₜ Q₂ Σ₂ Q₂ᵀ Pₜᵀ)
        primal_infeas = zero(R)
        @inbounds for (j, (nⱼ, rⱼ, r_pastⱼ, Wⱼ, Pⱼ, evⱼ)) in enumerate(zip(data.psds, data.r, r_past, data.W_psds, data.P_psds,
                                                                           cache.eigens))
            # note: we adjusted r such that it cannot exceed the side dimension of Xstar_psd, but we cannot do the same with
            # r_current and r_past, as only their sum has an upper bound.
            V = evⱼ[1]
            primal_infeas = min(primal_infeas, first(V.values))
            r_pastⱼ = min(r_pastⱼ, rⱼ)
            if iszero(r_pastⱼ)
                copyto!(Wⱼ, mastersolver.Wstar_psds[j])
                rmul!(Wⱼ, inv(tr(Wⱼ)))
                copyto!(Pⱼ, V.vectors)
            else
                γstarⱼ = max(mastersolver.γstars[j], zero(R)) # prevent numerical issues
                Sstareig = eigen!(mastersolver.Sstar_psds[j], @view(evⱼ[1].values[1:rⱼ]), evⱼ[5][:, 1:rⱼ], evⱼ[2], evⱼ[4])
                Q₁ = @view(Sstareig.vectors[:, end-r_pastⱼ+1:end]) # sorted in ascending order; we need the largest rₚ, but
                                                                   # the order doesn't really matter
                Q₂ = @view(Sstareig.vectors[:, 1:end-r_pastⱼ])
                Σ₂ = @view(Sstareig.values[1:end-r_pastⱼ])
                # Wⱼ = (γstar * Wⱼ + Pⱼ * Q₂ * Diagonal(Σ₂) * Q₂' * Pⱼ') / (γstar + tr(Σ₂))
                den = γstarⱼ + sum(v -> max(v, zero(R)), Σ₂) # also prevent numerical issues here
                #if den > sqrt(eps(R))
                    newpart = PackedMatrix(rⱼ, fill!(gettmp(cache, packedsize(rⱼ)), zero(R)), :L)
                    for (factor, newcol) in zip(Σ₂, eachcol(Q₂))
                        if factor > zero(R) # just to be sure
                            spr!(factor, newcol, newpart)
                        end
                    end
                    newpart_scaled = packed_scale!(newpart)
                    den = inv(den)
                    mul!(Wⱼ, cache.Pkrons[j], newpart_scaled, den, γstarⱼ * den)
                #end # else no update of W
                # Pⱼ = orth([V.vectors Pⱼ*Q₁])
                # for orthogonalization, we use QR to be numerically stable; unfortunately, this doesn't produce Q directly, so
                # we need another temporary. For consistency with the reference implementation, we put Pⱼ*Q₁ first (although it
                # uses orth, which is SVD-based).
                tmp = gettmp(cache, nⱼ, rⱼ)
                mul!(@view(tmp[:, 1:r_pastⱼ]), Pⱼ, Q₁)
                copyto!(@view(tmp[:, r_pastⱼ+1:end]), V.vectors)
                copyto!(Pⱼ, qr!(tmp).Q)
            end
        end
        # 15: if stopping criterion then
        #     16: quit
        isone(t) && continue

        # Our own dynamic strategy to increase ρ if necessary. If ρ was chosen too small, we will be able to reduce our primal
        # objective more and more at the expense of completely losing primal feasibility (because it cannot be achieved
        # anyway). This means that the constraint is (more than) active. Let's figure out the trace of the dual
        # Zⱼ = Cⱼ - ∑ᵢⱼ Aᵢⱼ yᵢ.
        if adaptiveρ
            trdual = zero(R)
            # Adapt the parameter ρ if necessary.
            for (aⱼ, Cⱼ) in zip(data.a_psds, data.C_psds)
                for i in PackedDiagonalIterator(Cⱼ)
                    trdual += Cⱼ[i] - dot(@view(aⱼ[:, i]), mastersolver.ystar)
                end
            end
            changed_ρ = trdual > R(1.1) * ρ
            if changed_ρ
                ρ *= R(2)
                specbm_adjust_penalty_subsolver!(cache.subsolver, ρ)
            end
        end

        # Iteration | Primal objective | Primal infeas | Dual infeas | Duality gap | Rel. accuracy | Rel. primal inf. | Rel. dual inf. | Rel. gap | Descent step | Consecutive null steps
        iszero(t % step) && @verbose_info(@sprintf("%9d | %16g | %13g | %11g | %11g | %13g | %16g | %14g | %11g | %12s | %22d",
            t, FΩ + offset, primal_infeas, dual_infeas, gap, rel_accuracy, rel_primal_infeas, rel_dual_infeas, rel_gap,
            has_descended, null_count), adaptiveρ ? @sprintf(" | %10g%s", trdual, changed_ρ ? " !" : "") : "")
        quality = max(rel_accuracy, rel_primal_infeas, rel_dual_infeas, rel_gap, -primal_infeas)
        if quality < ϵ
            status = :Optimal
            break
        end
        # 17: end if

        if !iszero(maxnodescent) && null_count ≥ maxnodescent
            status = :SlowProgress
            break
        end
    # 18: end for
    end

    specbm_finalize_primal_subsolver!(cache.subsolver)

    return SpecBMResult(status, FΩ + offset, data.Ω, mastersolver.ystar, t, quality, primal_infeas, dual_infeas, gap,
        rel_accuracy, rel_primal_infeas, rel_dual_infeas, rel_gap)
end

function specbm_setup_primal_subsolver end
function specbm_adjust_penalty_subsolver! end
function specbm_finalize_primal_subsolver! end
function specbm_primal_subsolve! end

if isdefined(Mosek, :appendafes)
    if VersionNumber(Mosek.getversion()) ≥ v"10.1.11"
        include("SpecBMMosek.jl")
    else
        @warn "The SpecBM method Mosek is not available: upgrade your Mosek distribution to at least version 10.1.11."
    end
end
include("SpecBMHypatia.jl")

if VERSION < v"1.10-"
    # identical to the implementation in SparseArrays, we just extend the allowed type for A, as this is already working
    # in Julia 1.10, the methods signatures were rewritten a lot and this is now supported natively.
    function LinearAlgebra.mul!(C::StridedVecOrMat, A::SparseArrays.SparseMatrixCSCView, B::SparseArrays.DenseInputVecOrMat,
        α::Number, β::Number)
        size(A, 2) == size(B, 1) || throw(DimensionMismatch())
        size(A, 1) == size(C, 1) || throw(DimensionMismatch())
        size(B, 2) == size(C, 2) || throw(DimensionMismatch())
        nzv = nonzeros(A)
        rv = rowvals(A)
        if β != 1
            β != 0 ? rmul!(C, β) : fill!(C, zero(eltype(C)))
        end
        for k in 1:size(C, 2)
            @inbounds for col in 1:size(A, 2)
                αxj = B[col,k] * α
                for j in nzrange(A, col)
                    C[rv[j], k] += nzv[j]*αxj
                end
            end
        end
        C
    end
end

@inline function direction_qp_primal_free!(mastersolver::SpecBMMastersolverData, data::SpecBMData, feasible::Bool, α::R,
    cache::SpecBMCache) where {R}
    invα = inv(α)
    # We need to (34): maximize dot(m, v) - dot(v, M, v) + const.
    #                      s.t. v = [γ; svec(S)]
    #                           γ ≥ 0, S ∈ 𝕊₊ʳ, γ + tr(S) ≤ ρ
    # Note that as we have multiple PSD blocks which we all treat separately (and not just as a single block-diagonal
    # constraint, we actually get multiple γ and multiple S matrices), though there is just one ρ.
    # Creating the data from the given parameters is detailed in C.1
    # We create a matrix Pkron (symmetrized Kronecked product) such that svec(Pᵀ W P) = Pkronᵀ*w, if w is the packed and scaled
    # vector of W. Note that due to the scaling, this is symmetric, so that svec(P U Pᵀ) = Pkron*u.
    # Pkronᵢ is packedsize(nᵢ) × packedsize(rᵢ)
    @inbounds @fastmath for (Pⱼ, Pkronⱼ) in zip(data.P_psds, cache.Pkrons)
        rows, cols = size(Pⱼ)
        colidx = 1
        for l in 1:cols
            rowidx = 1
            for k in 1:rows
                Pⱼkl = Pⱼ[k, l]
                Pkronⱼ[rowidx, colidx] = Pⱼkl^2
                rowidx += 1
                @simd for p in k+1:rows
                    Pkronⱼ[rowidx, colidx] = sqrt2 * Pⱼkl * Pⱼ[p, l]
                    rowidx += 1
                end
            end
            colidx += 1
            for q in l+1:cols
                rowidx = 1
                for k in 1:rows
                    Pⱼkl, Pⱼkq = Pⱼ[k, l], Pⱼ[k, q]
                    Pkronⱼ[rowidx, colidx] = sqrt2 * Pⱼkq * Pⱼkl
                    rowidx += 1
                    @simd for p in k+1:rows
                        Pkronⱼ[rowidx, colidx] = Pⱼkq * Pⱼ[p, l] + Pⱼkl * Pⱼ[p, q]
                        rowidx += 1
                    end
                end
                colidx += 1
            end
        end
    end
    # m₁ = q₁ - Q₁₃ Q₃₃⁻¹ q₃
    # q₁ = 2⟨Wⱼ, -α Ωⱼ + Cⱼ⟩
    # Q₃₁ = [⟨Wⱼ, Aᵢⱼ⟩]ᵢⱼ
    # q₃ = [2α(bᵢ - ⟨aᵢ, ω_free⟩ - ∑ⱼ ⟨Aᵢⱼ, Ωⱼ⟩) + 2(⟨c_free, aᵢ⟩ + ∑ⱼ ⟨Cⱼ, Aᵢⱼ⟩)
    # We can use Xstar_psd as temporaries for 2(-α Ωⱼ + Cⱼ)
    twoCminusαΩ = mastersolver.Xstar_psds
    mastersolver.xstar_psd .= R(2) .* (data.c_psd .- α .* data.ω_psd)
    cache.q₁ .= dot.(data.W_psds, twoCminusαΩ) # note that q₁ aliases m₁, so we already set the first part in m₁!
    mul!.(eachcol(cache.Q₃₁), data.a_psds, data.W_psds)
    if feasible
        copyto!(cache.q₃, cache.twoAc)
    else
        copyto!(cache.q₃, data.b)
        mul!(cache.q₃, data.A, data.Ω, R(-2) * α, R(2) * α)
        cache.q₃ .+= cache.twoAc
    end
    copyto!(mastersolver.ystar, cache.q₃) # we'll construct ystar successively, let's save q₃ for the moment
    ldiv!(cache.Q₃₃inv, cache.q₃) # now q₃ ← Q₃₃⁻¹ q₃
    mul!(cache.m₁, transpose(cache.Q₃₁), cache.q₃, -one(R), true)

    # m₂ = q₂ - Q₂₃ Q₃₃⁻¹ q₃
    # q₂ = (2vec(Pⱼᵀ (-α Ωⱼ + Cⱼ) Pⱼ))
    mul!.(cache.q₂s, transpose.(cache.Pkrons), twoCminusαΩ) # note that q₂s aliases m₂, so we already set the first part in m₂!
    # Q₃₂ = [vec(Pⱼᵀ Aᵢⱼ Pⱼ)ᵀ]ᵢⱼ
    mul!.(cache.Q₃₂s, data.a_psds, cache.Pkrons)
    mul!(cache.m₂, transpose(cache.Q₃₂), cache.q₃, -one(R), true) # q₃ already contains Q₃₃⁻¹ q₃

    # M₁₁ = Q₁₁ - Q₃₁ᵀ Q₃₃⁻¹ Q₃₁
    # Q₁₁ = Diag(⟨Wⱼ, Wⱼ⟩)
    tmpm = gettmp(cache, size(cache.Q₃₁)...)
    ldiv!(tmpm, cache.Q₃₃inv, cache.Q₃₁)
    mul!(cache.M₁₁, transpose(cache.Q₃₁), tmpm, -one(R), false)
    cache.Q₁₁ .+= LinearAlgebra.norm2.(data.W_psds) .^ 2 # Q₁₁ is a diagonal view into M₁₁

    # M₂₁ = Q₂₁ - Q₃₂ᵀ Q₃₃⁻¹ Q₃₁
    # Q₂₁ = Diag(svec(Pⱼᵀ Wⱼ Pⱼ)) - but this is a block diagonal for which there is no native support, so we use Vector{Vector}
    fill!(cache.M₂₁, zero(R))
    mul!.(cache.Q₂₁s, transpose.(cache.Pkrons), data.W_psds) # note that Q₂₁ aliases M₂₁, so we already set the first part!
    mul!(cache.M₂₁, transpose(cache.Q₃₂), tmpm, -one(R), true) # tmpm already contains the inverse part

    # M₂₂ = Q₂₂ - Q₃₂ᵀ Q₃₃⁻¹ Q₃₂
    # Q₂₂ = id_{Σr}
    tmpm = gettmp(cache, size(cache.Q₃₂)...)
    ldiv!(tmpm, cache.Q₃₃inv, cache.Q₃₂)
    mul!(cache.M₂₂, transpose(cache.Q₃₂), tmpm, -one(R), false)
    cache.Q₂₂ .+= one(R) # Q₂₂ is a diagonal view into M₂₂

    specbm_primal_subsolve!(mastersolver, cache)

    # Reconstruct y = Q₃₃⁻¹(q₃/2 - Q₃₁ γ - Q₃₂ svec(S))
    # Note that at this stage, we have already saved the original value of q₃ in y
    mul!(mastersolver.ystar, cache.Q₃₁, mastersolver.γstars, -one(R), inv(R(2)))
    mul!(mastersolver.ystar, cache.Q₃₂, mastersolver.sstar_psd, -one(R), true)
    ldiv!(cache.Q₃₃inv, mastersolver.ystar)
    # Reconstruct Wstarⱼ = γstarⱼ Wⱼ + Pⱼ Sstarⱼ Pⱼᵀ and Xstarⱼ = Ωⱼ + (Wstar - C + A*(ystar))/α
    copyto!(mastersolver.wstar_psd, data.w_psd)
    mul!.(mastersolver.Wstar_psds, cache.Pkrons, mastersolver.Sstar_psds, one(R), mastersolver.γstars)
    mastersolver.xstar_free .= .-data.c_free
    mastersolver.xstar_psd .= mastersolver.wstar_psd .- data.c_psd
    mul!(mastersolver.Xstar, data.At, mastersolver.ystar, invα, invα)
    # before we complete by adding Ω, calculate some feasibility quantifiers
    dfeasible_psd = (α * LinearAlgebra.norm2(mastersolver.xstar_psd))^2
    dfeasible_free = (α * norm(mastersolver.xstar_free))^2 # vector norm2 doesn't work on empty collections
    dfeasible = dfeasible_free + dfeasible_psd
    mastersolver.Xstar .+= data.Ω

    gap = abs(dot(data.b, mastersolver.ystar) - dot(data.c, mastersolver.Xstar))
    return dfeasible, dfeasible_free, dfeasible_psd, gap
end
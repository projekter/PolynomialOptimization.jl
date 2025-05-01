module FacialReduction

using MultivariatePolynomials, ..IntPolynomials, ..IntPolynomials.MultivariateExponents
using ..PolynomialOptimization: @assert, @verbose_info, FastKey, Problem, poly_problem
using ..Relaxation
using ..Solver: unique_outer_groupings

export facial_reduction!!

function M⁺(m::IntMonomialVector{Nr,Nc}) where {Nr,Nc}
    # M \ {(α + β)/2 : α, β ∈ M, α ≠ β}
    # We typically expect M⁺ to be rather small compared to m.
    indices = trues(length(m))
    α = Vector{Int}(undef, nvariables(eltype(m)))
    β = similar(α)
    γ = similar(α)
    @inbounds for (i, α) in enumerate(IntPolynomials.veciter(m, α))
        for β in Iterators.take(IntPolynomials.veciter(m, β), i -1)
            valid = true
            deg = 0
            for k in eachindex(α)
                γₖ = α[k] + β[k]
                if iszero(γₖ & true)
                    deg += (γ[k] = γₖ >> 1)
                else
                    valid = false
                    break
                end
            end
            valid || continue
            delEidx = IntPolynomials.MultivariateExponents.exponents_to_index(m.e, γ, deg)
            iszero(delEidx) && continue
            delidx = searchsortedlast(m, IntMonomial{Nr,Nc}(IntPolynomials.unsafe, m.e, delEidx, deg))
            if delidx > 0 && m[delidx].index == delEidx
                indices[delidx] = false
            end
        end
    end
    @inbounds return m[indices]
end

@inline function inΣhatofMperp(m::IntMonomial{Nr,Nc}, M::IntMonomialVector{Nr,Nc}, (mindegM, maxdegM)::Tuple{Int,Int},
                               newExponents::Vector{Int}) where {Nr,Nc}
    # p ∈ ̂Σ(M)⟂ = Σ(M)⟂, which is simply the set of all polynomials in ℝ_{2d} with exponents not in M + M
    # We need to figure out whether ∃m₁, m₂ ∈ M : exponents(m) = m₁ + m₂.
    # Clearly, this implies that degree(m₁) + mindegree(M) ≤ degree(m) ≤ degree(m₁) + maxdegree(M)
    # Therefore, we can restrict our search for m₁ to degree(m) - maxdegree(M) ≤ degree(m₁) ≤ degree(m) + mindegree(M)
    # and, more loosely, 2mindegree(M) ≤ degree(m) ≤ 2maxdegree(M)
    @assert(length(newExponents) == Nr + 2Nc)
    d = degree(m)
    drange = degree_range(M.e, d-maxdegM:d+mindegM)
    # It would be good to specify degree here; but we don't know whether the exponent set actually contains monomials of the
    # boundary degrees.
    subvector = @view(M[range(searchsortedfirst(M, IntMonomial{Nr,Nc}(IntPolynomials.unsafe, M.e, first(drange))),
                              searchsortedlast(M, IntMonomial{Nr,Nc}(IntPolynomials.unsafe, M.e, last(drange))))])
    @inbounds for m₁ in subvector
        failed = false
        degm₁ = degree(m₁)
        remdeg = d - degm₁
        for (j, (m₁ⱼ, mⱼ)) in enumerate(zip(exponents(m₁), exponents(m)))
            if (newExponents[j] = mⱼ - m₁ⱼ) < 0 || (remdeg -= newExponents[j]) < 0
                failed = true
                break
            end
        end
        failed && continue
        iszero(remdeg) || continue
        degm₂ = d - degm₁
        m₂idx = exponents_to_index(M.e, newExponents, degm₂)
        iszero(m₂idx) && continue
        insorted(IntMonomial{Nr,Nc}(IntPolynomials.unsafe, M.e, m₂idx, degm₂), subvector) && return false
    end
    return true
end

@inline function λindex(m::IntMonomial{Nr,Nc}, M⁺::IntMonomialVector{Nr,Nc}, newExponents::Vector{Int}) where {Nr,Nc}
    for (i, mᵢ) in enumerate(exponents(m))
        iszero(mᵢ & true) || return
        newExponents[i] = mᵢ >> 1
    end
    halfdeg = degree(m) >> 1
    M⁺Eidx = IntPolynomials.MultivariateExponents.exponents_to_index(M⁺.e, newExponents, halfdeg)
    iszero(M⁺Eidx) && return
    M⁺idx = searchsortedlast(M⁺, IntMonomial{Nr,Nc}(IntPolynomials.unsafe, M⁺.e, M⁺Eidx, halfdeg))
    M⁺idx > 0 && M⁺[M⁺idx].index == M⁺Eidx && return M⁺idx
    return
end

function λtoresult(λ::AbstractVector{<:Real}, M::AbstractVector{<:IntMonomialVector{Nr,Nc}},
                   M⁺::AbstractVector{<:IntMonomialVector{Nr,Nc}}, unusedλ::AbstractVector{Bool}) where {Nr,Nc}
    # The artificial constraint ∑λ = 1 is introduced just to exclude the case where everything is zero. But this means that we
    # could write any other positive constant on the r.h.s. Therefore, any λ that is never used anywhere in the problem may be
    # assigned an arbitrary strictly positive value, which allows to discard it - even if the solver assigned 0.
    λlast = 0
    changed = false
    result = similar(M)
    @inbounds for (j, (Mⱼ, M⁺ⱼ)) in enumerate(zip(M, M⁺))
        # TODO: write something analogous to keepat!!, using an iterator over deleted indices
        r = λlast+1:λlast+length(M⁺ⱼ)
        keeps = trues(length(Mⱼ))
        for (i, (λᵢ, unusedλᵢ)) in enumerate(zip(@view(λ[r]), @view(unusedλ[r])))
            if λᵢ > 1e-8 || unusedλᵢ
                keeps[searchsortedfirst(Mⱼ, M⁺ⱼ[i])] = false
                changed = true
            end
        end
        result[j] = IntPolynomials.keepat!!(Mⱼ, keeps)
        λlast += length(M⁺ⱼ)
    end
    return changed ? result : nothing
end

@doc raw"""
    facial_reduction(method, M, g; verbose=false, parameters...)

Performs facial reduction on a sums-of-squares problem. If just a single SOS condition is present, `M` may be an
`IntMonomialVector` and `g` a vector of `IntPolynomials`; for multiple conditions, `M` should be a vector of
`IntMonomialVector`s and `g` a matrix of `IntPolynomials` (with one column for each condition).
The vector(s) `M` contain initial candidates for the grouping while the polynomials `g` are fixed by the problem.

The form of the problem is
```math
\min_{\vec u \in \mathbb R^m} \vec \langle\vec w, \vec u\rangle \\
\text{such that} \\
g_{1, j}(\vec x) + \sum_{k = 2}^{m +1} u_k g_{k, j}(\vec x) \in \text{SOS}
```
and the algorithm follows [Permenter, Parillo (2014)](https://doi.org/10.1109/CDC.2014.7040427).
"""
function facial_reduction!!(method::Val, M::AbstractVector{<:IntMonomialVector{Nr,Nc,I}},
    gs::AbstractMatrix{<:IntPolynomial{<:Any,Nr,Nc,<:IntMonomialVector{Nr,Nc,I}}}; verbose::Bool=false, kwargs...) where {Nr,Nc,I<:Integer}
    length(M) == size(gs, 2) || throw(DimensionMismatch("Number of initial coordinate vectors and polynomials do not coincide"))
    result = M
    i = 0
    if verbose
        oldsize = sum(length, result, init=0)
    end
    while true
        @verbose_info("Facial reduction, iteration ", i += 1)
        frtime = @elapsed(new_result = facial_reduction!!(method, result, gs, M⁺.(result); verbose, kwargs...))
        if isnothing(new_result)
            @verbose_info("Found no further reduction in ", frtime, " seconds")
            return result
        else
            @verbose_info("Reduced ", oldsize, " elements to ", (oldsize = sum(length, new_result, init=0)), " in ",
                frtime, " seconds")
            result = new_result
        end
    end
end

facial_reduction!!(method::Val, M::IntMonomialVector{Nr,Nc,I}, g::AbstractVector{<:IntPolynomial{<:Any,Nr,Nc,I}};
    kwargs...) where {Nr,Nc,I<:Integer} = facial_reduction!(method, [M], [g;;]; kwargs...)

facial_reduction!!(method::Symbol, args...; kwargs...) = facial_reduction!!(Val(method), args...; kwargs...)
facial_reduction!!(M::AbstractVector, rest...; kwargs...) =
    facial_reduction!!(Val(default_reduction_method()), M, args...; kwargs...)

const reduction_methods = Symbol[]

function default_reduction_method()
    isempty(reduction_methods) &&
        error("No facial reduction method is available. Load a solver package that provides such a method (e.g., Mosek)")
    @inbounds return reduction_methods[begin]
end

include("./Constraints.jl")

end
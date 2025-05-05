module FacialReduction

using MultivariatePolynomials, ..FastVector, ..IntPolynomials, ..IntPolynomials.MultivariateExponents
using ..PolynomialOptimization: @assert, @capture, @inbounds, @verbose_info, FastKey, Problem, poly_problem
using ..Relaxation
using ..Solver: unique_outer_groupings, trisize

export facial_reduction

include("./Constraints.jl")

intsum(x...) = sum(x..., init=0)::Int

function M⁺(M::IntMonomialVector{Nr,0}, keeps::AbstractVector{Bool}, start::Integer, len::Integer) where {Nr}
    start > length(M) && return keeps # might happen due to roundoff errors
    α = Vector{Int}(undef, Nr)
    β = similar(α)
    γ = similar(α)
    @inbounds for (i, α) in zip(start:start+len-1, Iterators.drop(IntPolynomials.veciter(M, α), start -1))
        for β in Iterators.take(IntPolynomials.veciter(M, β), i -1)
            evendeg = true
            deg = 0
            for k in eachindex(α)
                γₖ = α[k] + β[k]
                if iszero(γₖ & true)
                    deg += (γ[k] = γₖ >> 1)
                else
                    evendeg = false
                    break
                end
            end
            evendeg || continue
            delEidx = IntPolynomials.MultivariateExponents.exponents_to_index(M.e, γ, deg)
            iszero(delEidx) && continue
            delidx = searchsortedlast(M, IntMonomial{Nr,0}(IntPolynomials.unsafe, M.e, delEidx, deg))
            if delidx > 0 && M[delidx].index == delEidx
                keeps[delidx] = false
            end
        end
    end
    return keeps
end

function M⁺(M::IntMonomialVector{Nr,0}; verbose::Bool=false) where {Nr}
    # M \ {(α + β)/2 : α, β ∈ M, α ≠ β}
    nthreads = Threads.nthreads()
    if (isone(nthreads) || length(M) < 128)
        @verbose_info("Using single-threaded setup, checking $(trisize(length(M))) combinations")
        @inbounds return M[M⁺(M, trues(length(M)), 1, length(M))]
    end
    # multithreading issues:
    # - different batch sizes - we parallelize the outer loop, but the inner has increasing length
    #   ↪ if a thread works from index start to start+δ-1, it covers δ*(2start + δ -3)/2 iterations
    #   → this is batchsize if δ = 3/2 - start + sqrt(2batchsize + (3/2 - start)^2)
    # - when we write to keeps, this must not be a bit array. At least UInt8 is necessary to avoid data races
    keeps = fill(true, length(M))
    totalsize = trisize(length(M) -1)
    nthreads = min(nthreads, max(1, totalsize >> 7)) # we don't want to have too small batchsizes, make 128 the minimum
    batchsize = div(totalsize, nthreads, RoundUp)
    @verbose_info("Using multi-threaded setup, each thread checking about $batchsize combinations")
    tasks = Vector{Task}(undef, nthreads)
    start = 1
    ccall(:jl_enter_threaded_region, Cvoid, ())
    try
        @inbounds for i in 1:nthreads
            δ = ceil(Int, 1.5 - start + sqrt(2batchsize + (1.5 - start)^2))
            tasks[i] = Threads.@spawn M⁺($M, $keeps, $start, $δ)
            start += δ
        end
        for task in tasks
            wait(task)
        end
    finally
        ccall(:jl_exit_threaded_region, Cvoid, ())
    end
    @inbounds return M[keeps]
end

function updateM⁺!(data::FacialReductionData; verbose::Bool=false)
    @verbose_info("Updating M⁺")
    upd_time = @elapsed(data.M⁺_obj .= M⁺.(data.M_obj; verbose))
    @verbose_info("├ objective (", upd_time, " seconds)")

    upd_time = @elapsed begin
        for (groupings, M⁺s) in zip(data.M_nonneg, data.M⁺_nonneg)
            M⁺s .= M⁺.(groupings; verbose)
        end
    end
    isempty(data.M_nonneg) || @verbose_info("├ nonnegative constraints (", upd_time, " seconds)")

    upd_time = @elapsed begin
        for (groupings, M⁺s) in zip(data.M_psd, data.M⁺_psd)
            M⁺s .= M⁺.(groupings; verbose)
        end
    end
    isempty(data.M_psd) || @verbose_info("├ PSD constraints (", conv_time, " seconds)")

    @verbose_info("└ End updating M⁺")
    return data
end

@inline function inΣhatofMperp(M::IntMonomialVector{Nr,0}, m::IntMonomial{Nr,0}, (mindegM, maxdegM)::Tuple{Int,Int},
                               newExponents::AbstractVector{Int}) where {Nr}
    # p ∈ ̂Σ(M)⟂ = Σ(M)⟂, which is simply the set of all polynomials in ℝ_{2d} with exponents not in M + M
    # We need to figure out whether ∃m₁, m₂ ∈ M : m = m₁ + m₂.
    # Clearly, this implies that degree(m₁) + mindegree(M) ≤ degree(m) ≤ degree(m₁) + maxdegree(M)
    # Therefore, we can restrict our search for m₁ to degree(m) - maxdegree(M) ≤ degree(m₁) ≤ degree(m) + mindegree(M)
    # and, more loosely, 2mindegree(M) ≤ degree(m) ≤ 2maxdegree(M)
    @assert(length(newExponents) ≥ Nr)
    d = degree(m)
    drange = degree_range(M.e, d-maxdegM:d+mindegM)
    # It would be good to specify degree here; but we don't know whether the exponent set actually contains monomials of the
    # boundary degrees.
    subvector = @view(M[range(searchsortedfirst(M, IntMonomial{Nr,0}(IntPolynomials.unsafe, M.e, first(drange))),
                              searchsortedlast(M, IntMonomial{Nr,0}(IntPolynomials.unsafe, M.e, last(drange))))])
    @inbounds for m₁ in subvector
        maybe_present = true
        degm₁ = degree(m₁)
        remdeg = d - degm₁
        for (j, (m₁ⱼ, mⱼ)) in enumerate(zip(exponents(m₁), exponents(m)))
            if (newExponents[j] = mⱼ - m₁ⱼ) < 0 || (remdeg -= newExponents[j]) < 0
                maybe_present = false
                break
            end
        end
        maybe_present || continue
        iszero(remdeg) || continue
        degm₂ = d - degm₁
        m₂idx = exponents_to_index(M.e, newExponents, degm₂)
        iszero(m₂idx) && continue
        insorted(IntMonomial{Nr,0}(IntPolynomials.unsafe, M.e, m₂idx, degm₂), subvector) && return false
    end
    return true
end

@inline function λindex(M⁺::IntMonomialVector{Nr,0}, m::IntMonomial{Nr,0}, newExponents::AbstractVector{Int}) where {Nr}
    @assert(length(newExponents) ≥ Nr)
    for (i, mᵢ) in enumerate(exponents(m))
        iszero(mᵢ & true) || return 0
        newExponents[i] = mᵢ >> 1
    end
    halfdeg = degree(m) >> 1
    M⁺Eidx = IntPolynomials.MultivariateExponents.exponents_to_index(M⁺.e, newExponents, halfdeg)
    iszero(M⁺Eidx) && return
    M⁺idx = searchsortedlast(M⁺, IntMonomial{Nr,0}(IntPolynomials.unsafe, M⁺.e, M⁺Eidx, halfdeg))
    @inbounds M⁺idx > 0 && M⁺[M⁺idx].index == M⁺Eidx && return M⁺idx
    return 0
end

function _truncate(M::IntMonomialVector{Nr,0}, M⁺::IntMonomialVector{Nr,0}, λ::AbstractVector{<:Real},
                   usedλ::AbstractVector{Bool}) where {Nr}
    # The artificial constraint ∑λ = 1 is introduced just to exclude the case where everything is zero. But this means that we
    # could write any other positive constant on the r.h.s. Therefore, any λ that is never used anywhere in the problem may be
    # assigned an arbitrary strictly positive value, which allows to discard it - even if the solver assigned 0.
    keeps = trues(length(M))
    changed = false
    for (i, (λᵢ, usedλᵢ)) in enumerate(zip(λ, usedλ))
        if λᵢ > 1e-8 || !usedλᵢ
            keeps[searchsortedfirst(M, M⁺[i])] = false
            changed = true
        end
    end
    return changed ? M[keeps] : M
end

function truncate!(data::FacialReductionData, λ::AbstractVector{<:Real}, usedλ::AbstractVector{Bool})
    @assert(length(λ) == length(usedλ))
    changed = false
    r = 0:0
    @inbounds for (i, (M, M⁺)) in enumerate(zip(data.M_obj, data.M⁺_obj))
        r = last(r)+1:last(r)+length(M⁺)
        newM = @views _truncate(M, M⁺, λ[r], usedλ[r])
        if newM !== M
            data.M_obj[i] = newM
            changed = true
        end
    end
    @inbounds for (groupings, M⁺s) in zip(data.M_nonneg, data.M⁺_nonneg)
        for (i, (M, M⁺)) in enumerate(zip(groupings, M⁺s))
            r = last(r)+1:last(r)+length(M⁺)
            newM = @views _truncate(M, M⁺, λ[r], usedλ[r])
            if newM !== M
                groupings[i] = newM
                changed = true
            end
        end
    end
    @inbounds for (groupings, M⁺s) in zip(data.M_psd, data.M⁺_psd)
        for (i, (M, M⁺)) in enumerate(zip(groupings, M⁺s))
            r = last(r)+1:last(r)+length(M⁺)
            newM = @views _truncate(M, M⁺, λ[r], usedλ[r])
            if newM !== M
                groupings[i] = newM
                changed = true
            end
        end
    end
    @assert(last(r) == length(λ))
    return changed
end

function facial_reduction! end

function _facial_reduction!(method::Val, data::FacialReductionData; verbose::Bool=false, kwargs...)
    i = 0
    @verbose_info("Starting facial reduction")
    while true
        if verbose
            println("Iteration #", i += 1)
            flush(stdout)
        end
        upd_time = @elapsed begin
            updateM⁺!(data; verbose)
            updateMonomials!(data; verbose)
        end
        @verbose_info("Iteration preprocessing done in ", upd_time, " seconds")
        frtime = @elapsed(changed = facial_reduction!(method, data; verbose, kwargs...))
        if !changed
            @verbose_info("Iteration time: ", frtime, " seconds - finished")
            return data
        else
            if verbose
                println("Iteration time: ", frtime, " seconds")
                show(stdout, "text/plain", data)
                flush(stdout)
            end
        end
    end
end

@doc raw"""
    facial_reduction!(method, relaxation::AbstractRelaxation; verbose=false, parameters...)

Performs facial reduction on a sums-of-squares problem.
The algorithm follows [Permenter, Parillo (2014)](https://doi.org/10.1109/CDC.2014.7040427).
"""
facial_reduction(@nospecialize(method::Val), relaxation::AbstractRelaxation; verbose::Bool=false, kwargs...) =
    # Relaxation.embed(
    #     @invoke(facial_reduction!(method::Val, FacialReductionData(relaxation; verbose)::FacialReductionData;
    #                               verbose, kwargs...)),
    #     Relaxation.groupings(relaxation);
    #     verbose
    # )
    Relaxation.embed(
        _facial_reduction!(method, FacialReductionData(relaxation; verbose); verbose, kwargs...),
        Relaxation.groupings(relaxation);
        verbose
    )

facial_reduction(method::Symbol, args...; kwargs...) = facial_reduction(Val(method), args...; kwargs...)
facial_reduction(relaxation::AbstractRelaxation; kwargs...) =
    facial_reduction(default_reduction_method(), relaxation; kwargs...)

const reduction_methods = Symbol[]

function default_reduction_method()
    isempty(reduction_methods) &&
        error("No facial reduction method is available. Load a solver package that provides such a method (e.g., Mosek)")
    @inbounds return reduction_methods[begin]
end

end
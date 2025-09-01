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
    isempty(data.M_psd) || @verbose_info("├ PSD constraints (", upd_time, " seconds)")

    @verbose_info("└ End updating M⁺")
    return data
end

@inline inΣhatofMperp(M²::IntMonomialVector{Nr,0}, m::IntMonomial{Nr,0}) where {Nr} = !(m ∈ M²)
    # p ∈ ̂Σ(M)⟂ = Σ(M)⟂, which is simply the set of all polynomials in ℝ_{2d} with exponents not in M + M, i.e. the M² set.

@inline function λindex(M⁺::IntMonomialVector{Nr,0}, m::IntMonomial{Nr,0}, newExponents::AbstractArray{Int}) where {Nr}
    @assert(length(newExponents) ≥ Nr)
    iszero(degree(m) & true) || return 0
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
    changed || return M
    result = M[keeps]
    return result
end

function truncate!(data::FacialReductionData, λ::AbstractVector{<:Real}, usedλ::AbstractVector{Bool})
    # We want to get reproducible behavior, so an interruption after an iteration was displayed should always yield the same
    # result. Hence, we disable interrupts during truncation.
    @assert(length(λ) == length(usedλ))
    changed = false
    Base.sigatomic_begin()
    try
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
    finally
        Base.sigatomic_end()
    end
    return changed
end

function facial_reduction! end

function facial_reduction!(method::Val, data::FacialReductionData; verbose::Bool=false, maxiter::Integer=0, kwargs...)
    i = 1
    @verbose_info("Starting facial reduction")
    try
        while true
            @verbose_info("Iteration #", i)
            upd_time = @elapsed begin
                updateM⁺!(data; verbose)
                updateM²!(data; verbose)
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
            if (i += 1) == maxiter
                @verbose_info("Terminating due to maximum number of iterations")
                return data
            end
        end
    catch e
        if e isa InterruptException
            @verbose_info("Interrupted facial reduction, returning data of last iteration")
            return data
        end
        rethrow(e)
    end
end

@doc raw"""
    facial_reduction!(method, relaxation::AbstractRelaxation; verbose=false, [maxiter], parameters...)

Performs facial reduction on a sums-of-squares problem.
The algorithm follows [Permenter, Parillo (2014)](https://doi.org/10.1109/CDC.2014.7040427); it is iterative and will terminate
if no further reduction could be found in one iteration or if `maxiter` was reached. Cancelling the algorithm will return the
data from the last complete iteration.
"""
facial_reduction(@nospecialize(method::Val), relaxation::AbstractRelaxation; verbose::Bool=false, kwargs...) =
    Relaxation.embed!(
        @invoke(facial_reduction!(method::Val, FacialReductionData(relaxation; verbose)::FacialReductionData;
                                  verbose, kwargs...)),
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
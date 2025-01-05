export ExponentsAll

"""
    ExponentsAll{N,I}()

Represents an index range for a fixed number of variables `N` without any degree bound.
"""
struct ExponentsAll{N,I<:Integer} <: AbstractExponentsUnbounded{N,I} end

# Use a trick. ExponentsAll should be a singleton, but we also need it to contain the cache of the exponents, and it must
# be fast - the address of the cache should be inlined. We cannot make it a field, for then it would no longer be a singleton.
# So instead, we have one generated function that, per singleton, creates the cache upon compile-time as a Ref and always
# returns this value statically. This appears to work despite the restrictions on generated functions and heap allocations.

@inline @generated function index_counts(::Type{Ref}, ::ExponentsAll{N,I}) where {N,I<:Integer}
    N isa Integer || throw(MethodError(index_counts{N,I}, ()))
    N < 0 && error("The number of variables must be nonnegative")
    counts = Ref{Matrix{I}}(ones(I, 1, N +1)) # create an initial matrix in the beginning so that it is always guaranteed to
                                              # have N+1 columns (for convenience, fill the last column with 1, which
                                              # corresponds to zero variables) and at least one row
    # TODO: In Julia 1.11, arrays just wrap the memory. We might be able to exploit this and avoid the extra reference.
    :(return $counts)
end

@inline index_counts(::Unsafe, e::ExponentsAll) = index_counts(Ref, e)[]
@inline _has_index_counts(e::ExponentsAll, degree::Integer) = size(index_counts(Ref, e)[], 1) > degree
@generated function _calc_index_counts!(e::ExponentsAll{N,I}, degree::Integer) where {N,I<:Integer}
    # put this check in the generation part, but not outside of the function - we don't want it to be fixed during
    # precompilation!
    multithread = Threads.nthreads() > 1
    if multithread
        lock = Threads.Condition()
    end
    quote
        $((multithread ? (:(lock($lock)),) : ())...)
        counts = index_counts(Ref, e)
        oldcounts = counts[]
        oldrows = size(oldcounts, 1)
        if $(!multithread) || oldrows ≤ degree # maybe we already increased before the lock took place
            newcounts = Matrix{I}(undef, degree +1, N +1)
            @inbounds fill!(@view(newcounts[:, N+1]), one(I))
            for j in N:-1:1
                @inbounds copyto!(@view(newcounts[:, j]), @view(oldcounts[:, j]))
                for i in oldrows+1:degree+1
                    # Explicitly, we'd have
                    # newcol[i] = binomial(I(i -1) + I(nvars), I(nvars))
                    # But this is wasteful, we just can exploit Pascal's triangle for the calculation:
                    # binomial(i -1 + nvars, nvars) = binomial(i -2 + nvars, nvars -1) + binomial(i -2 + nvars, nvars)
                    @inbounds newcounts[i, j] = newcounts[i, j+1] + newcounts[i-1, j]
                    # as the last column and first row is pre-filled, this will always work
                end
            end
            counts[] = newcounts
        end
        $((multithread ? (:(unlock($lock)),) : ())...)
        return
    end
end

function _exponents_to_index(e::ExponentsAll{N,I}, exponents, degree::Int, report_lastexp) where {N,I<:Integer}
    iszero(degree) && return isnothing(report_lastexp) ? one(I) : (one(I), 0)
    counts, success = index_counts(e, degree)
    @assert(success)
    index::I = @inbounds counts[degree+1, 1]
    lastexp = -1
    @inbounds for (i, vardeg) in zip(2:(isnothing(report_lastexp) ? N : report_lastexp+1), exponents)
        lastexp = vardeg
        i == N +1 && break # just for report_lastexp, where we must visit the last exponent also.

        # We still need to get mondeg for the total degree, but the current variable only has vardeg. Skip over all the
        # exponents where the current variable had a higher degree - these are given by the total number of exponents where the
        # variables to the right of the current one have degree exactly mondeg-(vardeg+1), mondeg-(vardeg+2), ....
        degree > vardeg && (index -= counts[degree-vardeg, i])
        if iszero(degree -= vardeg) # shortcut
            if !isnothing(report_lastexp) && i ≤ report_lastexp
                lastexp = 0
            end
            break
        end
    end
    return isnothing(report_lastexp) ? index : (index, lastexp)
end

@inline function degree_from_index(::Unsafe, e::ExponentsAll{<:Any,I}, index::I) where {I<:Integer}
    @inbounds counts = @view(index_counts(unsafe, e)[:, 1])
    return searchsortedfirst(counts, index) -1
end

function degree_from_index(e::ExponentsAll{N,I}, index::I) where {N,I<:Integer}
    @inbounds counts = @view(index_counts(unsafe, e)[:, 1])
    @inbounds if counts[end] < index
        # Let's calculate the indices manually with binomial. In this way, we can fairly cheaply get a good guess, whereas
        # extending the cache would mean a lot of potentially superfluous allocations.
        lo = length(counts) -1
        hi = 2max(1, length(counts) -1)
        while binomial(I(hi) + I(N), I(N)) < index
            lo = hi
            hi *= 2
        end
        while lo < hi
            mid = Base.midpoint(lo, hi)
            midindex = binomial(I(mid) + I(N), I(N))
            if midindex < index
                lo = mid +1
            elseif midindex > index
                hi = mid
            else
                index_counts(e, mid)
                return mid
            end
        end
        # now populate the cache - we need one more to have an upper bound
        index_counts(e, lo +1) # TODO: check
        return lo
    else
        return searchsortedfirst(counts, index) -1
    end
end

function exponents_from_index(e::ExponentsAll{<:Any,I}, index::I, degree::Int) where {I<:Integer}
    index > zero(I) || throw(BoundsError(e, index))
    counts, success = index_counts(e, degree) # initialize the cache
    @assert(success)
    @inbounds if (iszero(degree) && index > 1) ||
        (!iszero(degree) && (counts[degree, 1] ≥ index || counts[degree+1, 1] < index))
        throw(ArgumentError("Index $index does not have degree $degree"))
    end
    return ExponentIndices(e, index, degree)
end

function Base.iterate(efi::ExponentIndices{I,<:ExponentsAll{<:Any,I}}) where {I<:Integer}
    counts = index_counts(unsafe, efi.parent)
    degree = efi.degree
    @inbounds return iterate(efi, (degree, 2, iszero(degree) ? efi.index : efi.index - counts[degree, 1]))
end

function Base.iterate(efi::ExponentIndices{I,ExponentsAll{N,I}}, (degree, i, index)::Tuple{Int,Int,I}) where {N,I<:Integer}
    counts = index_counts(unsafe, efi.parent) # must not be passed in the state - then Julia would have to do an allocation for
                                              # the return type
    @inbounds if i ≤ N
        iszero(degree) && return 0, (degree, i +1, index)
        # Obtain the degree for index i-1 by looking at how large the subspaces to the right of i-1 are when we fix i-1 to
        # a certain value - we are in the last subspace that still fits. The sizes of the exact subspaces are given by the
        # differences of adjacent rows in counts (as counts is accumulating). We must start from the smallest degree for
        # i-1, which means the largest degree for the subspace to the right. So we must accumulate in the opposite order
        # than what is stored in counts. But using the telescopic sum, these accumulations simplify.
        tmp = counts[degree+1, i]
        remainingdeg = searchsortedlast(@view(counts[:, i]), tmp - index)
        degᵢ₋₁ = degree - remainingdeg
        index -= tmp - counts[remainingdeg+1, i]
        return degᵢ₋₁, (remainingdeg, i +1, index)
    elseif i == N +1
        return degree, (0, i +1, zero(index))
    else
        return nothing
    end
end

function iterate!(::Unsafe, v::AbstractVector{Int}, e::ExponentsAll)
    @inbounds begin
        while true
            # This is not a loop at all, we only go through it once, but we need to be able to leave the block at multiple
            # positions. If we do it with @goto, as would be proper, Julia begins to box all our arrays.

            # find the next exponent that can be decreased
            i = findlast(>(0), v)
            isnothing(i) && break

            # we must increment the exponents to the left by 1 in total
            isone(i) && break
            v[i-1] += 1
            # this implies that we reset everything to the right of the increment to its minimum and then compensate for all
            # the reductions by increasing the exponents again
            δ = v[i] -1
            v[i] = 0
            v[end] += δ
            return true
        end
        # there's still hope: we can perhaps go to the next degree
        deg = sum(v, init=0) +1
        fill!(@view(v[1:end-1]), 0)
        @inbounds v[end] = deg
        return true
    end
end
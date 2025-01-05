export ExponentsDegree

mutable struct ExponentsDegree{N,I<:Integer} <: AbstractExponentsDegreeBounded{N,I}
    const mindeg::Int
    const maxdeg::Int
    counts::Matrix{I}

    @doc """
        ExponentsDegree{N,I}(mindeg::Integer, maxdeg::Integer)

    Represents an exponent range that is restricted by a bound on the total degree.
    """
    function ExponentsDegree{N,I}(mindeg::Integer, maxdeg::Integer) where {N,I<:Integer}
        0 ≤ N || throw(MethodError(ExponentsDegree{N}, (mindeg, maxdeg)))
        0 ≤ mindeg ≤ maxdeg || throw(ArgumentError("Invalid degree specification"))
        new{N,I}(mindeg, maxdeg)
    end
end

ExponentsDegree{N,I}(range::AbstractUnitRange) where {N,I<:Integer} = ExponentsDegree{N,I}(first(range), last(range))

Base.:(==)(e1::E, e2::E) where {E<:ExponentsDegree} = e1.mindeg == e2.mindeg && e1.maxdeg == e2.maxdeg

function _calc_index_counts!(E::ExponentsDegree{N,I}) where {N,I<:Integer}
    # This is not exactly correct if we think about the very last column (corresponding to zero variables), which should be
    # zero if a nonzero mindeg is set. But this last column is not used anyway (apart from the construction of the matrix), so
    # instead of needlessly duplicating data, we just use a reference to the current unbounded one. This is very efficient -
    # unless the unbounded one then grows, because we'll keep the reference alive unneccessarily.
    E.counts, success = index_counts(ExponentsAll{N,I}(), E.maxdeg)
    @assert(success)
    return
end

function _exponents_to_index(e::ExponentsDegree{N,I}, exponents, degree::Int, report_lastexp) where {N,I<:Integer}
    e.mindeg ≤ degree ≤ e.maxdeg || return isnothing(report_lastexp) ? zero(I) : (zero(I), -1)
    iszero(degree) && return isnothing(report_lastexp) ? one(I) : (one(I), 0)
    counts, success = @inbounds index_counts(e, degree)
    @assert(success)
    # Our index starts with the last exponent that has a degree ≤ the required
    index::I = @inbounds counts[degree+1, 1]
    iszero(e.mindeg) || (index -= @inbounds counts[e.mindeg, 1])
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

function Base.iterate(efi::ExponentIndices{I,ExponentsDegree{N,I}}, (degree, i, index)::Tuple{Int,Int,I}) where {N,I<:Integer}
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

function iterate!(::Unsafe, v::AbstractVector{Int}, e::ExponentsDegree)
    # we must assume that a valid vector is provided, i.e. sum(v) ≥ e.mindeg
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
        deg > e.maxdeg && return false
        fill!(@view(v[1:end-1]), 0)
        @inbounds v[end] = deg
        return true
    end
end
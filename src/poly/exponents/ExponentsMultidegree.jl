export ExponentsMultideg

mutable struct ExponentsMultideg{N,I<:Integer,Vmin,Vmax} <: AbstractExponentsDegreeBounded{N,I}
    const mindeg::Int
    const maxdeg::Int
    const minmultideg::Vmin
    const maxmultideg::Vmax
    const Σminmultideg::Int
    const Σmaxmultideg::Int
    counts::Matrix{I}

    @doc """
        ExponentsMultideg{N,I}(mindeg::Integer, maxdeg::Integer,
            minmultideg::AbstractVector, maxmultideg::AbstractVector)

    Represents an exponent range that is restricted both by a global bound on the degree and by individual bounds on the
    variable degrees. Note that the vectors must not be used afterwards, and the constructor may clip maxmultideg to be no
    larger than maxdeg in each entry.
    """
    function ExponentsMultideg{N,I}(mindeg::Integer, maxdeg::Integer, minmultideg::Vmin, maxmultideg::Vmax) where
        {N,I<:Integer,Vmin<:AbstractVector{<:Integer},Vmax<:AbstractVector{<:Integer}}
        0 < N || throw(MethodError(ExponentsDegree{N}, (mindeg, maxdeg)))
        0 ≤ mindeg ≤ maxdeg || throw(ArgumentError("Invalid degree specification"))
        length(minmultideg) == length(maxmultideg) == N || throw(ArgumentError("Unsuitable multidegree lengths"))
        Σminmultideg = 0
        Σmaxmultideg = 0
        @inbounds for i in eachindex(minmultideg, maxmultideg)
            maxmultideg[i] > maxdeg && (maxmultideg[i] = maxdeg)
            mindeg < minmultideg[i] && (mindeg = minmultideg[i])
            0 ≤ minmultideg[i] ≤ maxmultideg[i] || throw(ArgumentError("Invalid multidegree specification"))
            Σmaxmultideg = Base.Checked.checked_add(Σmaxmultideg, Int(maxmultideg[i]))
            Σminmultideg += Int(minmultideg[i])
        end
        mindeg = max(mindeg, Σminmultideg)
        maxdeg = min(maxdeg, Σmaxmultideg)
        0 ≤ mindeg ≤ maxdeg || throw(ArgumentError("Invalid multidegree specification"))
        new{N,I,Vmin,Vmax}(mindeg, maxdeg, minmultideg, maxmultideg, Σminmultideg, Σmaxmultideg)
    end
end

ExponentsMultideg{N,I}(range::AbstractUnitRange, minmultideg::Vmin, maxmultideg::Vmax) where {N,I<:Integer,Vmin<:AbstractVector{<:Integer},Vmax<:AbstractVector{<:Integer}} =
    ExponentsMultideg{N,I}(first(range), last(range), minmultideg, maxmultideg)

Base.:(==)(e1::ExponentsMultideg{N,I}, e2::ExponentsMultideg{N,I}) where {N,I<:Integer} =
    e1.mindeg == e2.mindeg && e1.maxdeg == e2.maxdeg && e1.minmultideg == e2.minmultideg && e1.maxmultideg == e2.maxmultideg

function _calc_index_counts!(e::ExponentsMultideg{N,I}) where {N,I<:Integer}
    maxdeg = e.maxdeg
    if e.mindeg > e.maxdeg
        e.counts = Matrix{I}(undef, 0, N +1)
        return
    end
    minmultideg = e.minmultideg
    maxmultideg = e.maxmultideg
    counts = zeros(I, maxdeg +1, N +1)
    @inbounds @views begin
        if iszero(e.mindeg)
            fill!(counts[:, N+1], 1)
        end
        nextround = counts[:, N]
        for (occₖ, indₖ) in enumerate(minmultideg[N]+1:maxmultideg[N]+1)
           nextround[indₖ] = occₖ
        end
        fill!(nextround[maxmultideg[N]+2:maxdeg+1], maxmultideg[N]-minmultideg[N]+1)
        for j in N-1:-1:1
            lastround = nextround
            nextround = counts[:, j]
            for degᵢ in minmultideg[j]:maxmultideg[j]
                for (degₖ, occₖ) in zip(Iterators.countfrom(0), lastround)
                    newdeg = degᵢ + degₖ
                    newdeg > maxdeg && break
                    nextround[newdeg+1] += occₖ
                end
            end
        end
    end
    e.counts = counts
    return
end

function _exponents_to_index(e::ExponentsMultideg{N,I}, exponents, degree::Int, report_lastexp) where {N,I<:Integer}
    e.mindeg ≤ degree ≤ e.maxdeg || return zero(I)
    iszero(degree) && return isnothing(report_lastexp) ? one(I) : (one(I), 0)
    counts, success = @inbounds index_counts(e, degree)
    @assert(success)
    # Our index starts with the last exponent that has a degree ≤ the required
    index::I = @inbounds counts[degree+1, 1]
    iszero(e.mindeg) || (index -= @inbounds counts[e.mindeg, 1])
    Σminmultideg_rem = e.Σminmultideg
    Σmaxmultideg_rem = e.Σmaxmultideg
    lastexp = -1
    @inbounds for (i, vardeg, minmultideg, maxmultideg) in zip(2:(isnothing(report_lastexp) ? N : report_lastexp+1), exponents,
                                                               e.minmultideg, e.maxmultideg)
        minmultideg ≤ vardeg ≤ maxmultideg || return isnothing(report_lastexp) ? zero(I) : (zero(I), vardeg)
        lastexp = vardeg
        i == N +1 && break # just for report_lastexp, where we must visit the last exponent also.

        # We still need to get mondeg for the total degree, but the current variable only has vardeg. Skip over all the
        # exponents where the current variable had a higher degree - these are given by the total number of exponents where the
        # variables to the right of the current one have degree exactly mondeg-(vardeg+1), mondeg-(vardeg+2), ...,
        # mondeg-maxmultideg.
        degree > vardeg && (index -= counts[degree-vardeg, i])
        degree > maxmultideg && (index += counts[degree-maxmultideg, i])
        degree -= vardeg
        Σminmultideg_rem -= minmultideg
        Σmaxmultideg_rem -= maxmultideg
    end
    # We take the slightly more expensive approach to accumulate the remaining multidegs instead of just accessing
    # last(multideg). This allows exponents to have less entries than required, calculating the last index in this subspace.
    # The check is redundant if exponents was complete.
    Σminmultideg_rem ≤ degree ≤ Σmaxmultideg_rem || return isnothing(report_lastexp) ? zero(I) : (zero(I), lastexp)
    return isnothing(report_lastexp) ? index : (index, lastexp)
end

function Base.iterate(efi::ExponentIndices{I,<:ExponentsMultideg{N,I}}, (degree, i, index)::Tuple{Int,Int,I}) where {N,I<:Integer}
    counts = index_counts(unsafe, efi.parent) # must not be passed in the state - then Julia would have to do an allocation for
                                              # the return type
    minmultideg, maxmultideg = efi.parent.minmultideg, efi.parent.maxmultideg
    @inbounds if i ≤ N
        iszero(degree) && return 0, (degree, i +1, index)
        # Obtain the degree for index i-1 by looking at how large the subspaces to the right of i-1 are when we fix i-1 to
        # a certain value - we are in the last subspace that still fits. The sizes of the exact subspaces are given by the
        # differences of adjacent rows in counts (as counts is accumulating). We must start from the smallest degree for
        # i-1, which means the largest degree for the subspace to the right. So we must accumulate in the opposite order
        # than what is stored in counts. But using the telescopic sum, these accumulations simplify.
        lowerbound = max(1, degree - maxmultideg[i-1])
        tmp = counts[degree-minmultideg[i-1]+1, i]
        remainingdeg = searchsortedlast(@view(counts[lowerbound:degree-minmultideg[i-1], i]), tmp - index) + lowerbound -1
        degᵢ₋₁ = degree - remainingdeg
        index -= tmp - counts[remainingdeg+1, i]
        return degᵢ₋₁, (remainingdeg, i +1, index)
    elseif i == N +1
        return degree, (0, i +1, zero(index))
    else
        return nothing
    end
end

function iterate!(::Unsafe, v::AbstractVector{Int}, e::ExponentsMultideg{N}) where {N}
    @inbounds begin
        minmultideg, maxmultideg = e.minmultideg, e.maxmultideg
        while true
            # This is not a loop at all, we only go through it once, but we need to be able to leave the block at multiple
            # positions. If we do it with @goto, as would be proper, Julia begins to box all our arrays.

            # find the next exponent that can be decreased
            found = false
            local i
            for outer i in lastindex(v):-1:firstindex(v)
                if v[i] > minmultideg[i]
                    found = true
                    break
                end
            end
            found || break
            # we must increment the exponents to the left by 1 in total
            found = false
            local j
            for outer j in i-1:-1:firstindex(v)
                if v[j] < maxmultideg[j]
                    found = true
                    break
                end
            end
            found || break

            v[j] += 1
            # this implies that we reset everything to the right of the increment to its minimum and then compensate for all
            # the reductions by increasing the exponents again
            δ = -1
            for k in j+1:i
                δ += v[k] - minmultideg[k]
            end
            copyto!(v, j +1, minmultideg, j +1, i - j)
            exponents_increment_right!(v, e, δ, j +1) && return true
        end
        # there's still hope: we can perhaps go to the next degree
        deg = sum(v, init=0) +1
        deg > e.maxdeg && return false
        copyto!(v, minmultideg)
        return exponents_increment_right!(v, e, deg - e.Σminmultideg, firstindex(v))
    end
end

function exponents_increment_right!(v::AbstractVector{Int}, e::ExponentsMultideg{N}, δ, from) where {N}
    @assert(δ ≥ 0 && from ≥ 0)
    maxmultideg = e.maxmultideg
    i = N
    @inbounds while δ > 0 && i ≥ from
        δᵢ = maxmultideg[i] - v[i]
        if δᵢ ≥ δ
            v[i] += δ
            return true
        else
            v[i] = maxmultideg[i]
            δ -= δᵢ
        end
        i -= 1
    end
    return iszero(δ)
end
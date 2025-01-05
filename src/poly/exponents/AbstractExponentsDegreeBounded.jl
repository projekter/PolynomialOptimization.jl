export AbstractExponentsDegreeBounded

"""
    AbstractExponentsDegreeBounded{N,I} <: AbstractExponents{N,I}

Abstract supertype for finite collections of multivariate exponents, bounded by their degrees. These collections have a length;
they also provide at least the fields `mindeg` and `maxdeg` that describe (a superset of) the covered degree range. Their
cache is always initialized completely when required.
"""
abstract type AbstractExponentsDegreeBounded{N,I<:Integer} <: AbstractExponents{N,I} end
# AbstractExponentsDegreeBounded descendants must implement
# - _has_index_counts(::AbstractExponents) that checks whether the cache matrix is fully initialized (should be fast)
# - _calc_index_counts!(::AbstractExponents) that initializes the full cache matrix.

Base.IteratorSize(::Type{<:AbstractExponentsDegreeBounded}) = Base.HasLength()
"""
    length(unsafe, e::AbstractExponentsDegreeBounded)

Unsafe variant of [`length`](@ref length(::AbstractExponentsDegreeBounded)): requires the cache to be set up
correctly, else the behavior is undefined.
"""
function Base.length(::Unsafe, e::AbstractExponentsDegreeBounded)
    counts = index_counts(unsafe, e)
    @inbounds return counts[e.maxdeg+1, 1] - (iszero(e.mindeg) ? zero(eltype(counts)) : counts[e.mindeg, 1])
end
"""
    length(e::AbstractExponentsDegreeBounded)

Returns the total number of exponents present in `e`.
"""
function Base.length(e::AbstractExponentsDegreeBounded)
    counts, success = index_counts(e, e.maxdeg)
    @assert(success)
    @inbounds return counts[e.maxdeg+1, 1] - (iszero(e.mindeg) ? zero(eltype(counts)) : counts[e.mindeg, 1])
end

Base.lastindex(e::AbstractExponentsDegreeBounded{<:Any,I}) where {I<:Integer} = I(length(e))

Base.:(==)(e1::AbstractExponentsDegreeBounded{N,I}, e2::AbstractExponentsDegreeBounded{N,I}) where {N,I<:Integer} =
    index_counts(e1, e1.maxdeg) == index_counts(e2, e2.maxdeg)

function Base.iterate(e::AbstractExponentsDegreeBounded{<:Any,I}) where {I<:Integer}
    counts, success = index_counts(e, e.mindeg) # ExponentIndices requires the cache to be set up
    iszero(e.mindeg) && return ExponentIndices(e, one(I), e.mindeg), (one(I), e.mindeg, one(I))
    @assert(success)
    @inbounds return ExponentIndices(e, one(I), e.mindeg), (one(I), e.mindeg, counts[e.mindeg+1, 1] - counts[e.mindeg, 1])
end
function Base.iterate(e::AbstractExponentsDegreeBounded{<:Any, I}, (index, degree, remainingdeg)::Tuple{I,Int,I}) where {I<:Integer}
    if isone(remainingdeg)
        degree == e.maxdeg && return nothing
        counts, success = index_counts(e, degree +1)
        @assert(success)
        @inbounds return ExponentIndices(e, index + one(I), degree +1),
            (index + one(I), degree +1, counts[degree+2, 1] - counts[degree+1, 1])
    else
        return ExponentIndices(e, index + one(I), degree), (index + one(I), degree, remainingdeg - one(I))
    end
end

Base.IteratorSize(::Type{<:ExponentsVectorIterator{<:Any,<:AbstractExponentsDegreeBounded}}) = Base.HasLength()
Base.length(ei::ExponentsVectorIterator{<:Any,<:AbstractExponentsDegreeBounded}) = length(unsafe, ei.e)

Base.@assume_effects :total !:inaccessiblememonly index_counts(::Unsafe, e::AbstractExponentsDegreeBounded) =
    (@inline; e.counts)

Base.@assume_effects :nothrow function index_counts(e::AbstractExponentsDegreeBounded, degree::Integer)
    @inline
    @inline(_has_index_counts(e)) || @noinline(_calc_index_counts!(e))
    return e.counts, size(e.counts, 1) > degree
end

@inline _has_index_counts(e::AbstractExponentsDegreeBounded) = isdefined(e, :counts)

@inline function degree_from_index(::Unsafe, e::AbstractExponentsDegreeBounded{<:Any,I}, index::I) where {I<:Integer}
    @inbounds counts = @view(index_counts(unsafe, e)[:, 1])
    iszero(e.mindeg) || (index += @inbounds counts[e.mindeg, 1])
    return searchsortedfirst(counts, index) -1
end

Base.@propagate_inbounds function degree_from_index(e::AbstractExponentsDegreeBounded{<:Any,I}, index::I) where {I<:Integer}
    counts = let ic=index_counts(e, e.maxdeg)
        @assert(ic[2])
        @inbounds @view(ic[1][:, 1])
    end
    iszero(e.mindeg) || (index += @inbounds counts[e.mindeg, 1])
    return searchsortedfirst(counts, index) -1
end

function exponents_from_index(e::AbstractExponentsDegreeBounded{<:Any,I}, index::I, degree::Int) where {I<:Integer}
    index > zero(I) || throw(BoundsError(e, index))
    e.mindeg ≤ degree ≤ e.maxdeg || throw(BoundsError(e, index))
    counts, success = index_counts(e, degree) # initialize the cache
    @assert(success)
    allindex = index
    iszero(e.mindeg) || (allindex += @inbounds counts[e.mindeg, 1])
    @inbounds if (iszero(degree) && allindex > 1) ||
        (!iszero(degree) && (counts[degree, 1] ≥ allindex || counts[degree+1, 1] < allindex))
        throw(ArgumentError("Index $index does not have degree $degree"))
    end
    return ExponentIndices(e, index, degree)
end

function Base.iterate(efi::ExponentIndices{I,<:AbstractExponentsDegreeBounded{<:Any,I}}) where {I<:Integer}
    parent = efi.parent
    counts = index_counts(unsafe, parent)
    degree = efi.degree
    @inbounds return iterate(efi, (degree, 2, (iszero(degree) ? efi.index : efi.index - counts[degree, 1]) +
                                              (iszero(parent.mindeg) ? zero(I) : counts[parent.mindeg, 1])))
end

include("./ExponentsDegree.jl")
include("./ExponentsMultidegree.jl")
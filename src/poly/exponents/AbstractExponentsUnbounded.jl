export AbstractExponentsUnbounded

"""
    AbstractExponentsUnbounded{N,I} <: AbstractExponents{N,I}

Abstract supertype for unbounded collections of multivariate exponents. These collections do not have a length; they are
infinite. Their cache is always initialized incrementally as required.
"""
abstract type AbstractExponentsUnbounded{N,I<:Integer} <: AbstractExponents{N,I} end
# AbstractExponentsUnbounded descendants must implement
# - _has_index_counts(::AbstractExponentsUnbounded, degree::Integer) that checks whether the cache matrix is valid for the
#   given degree (should be fast)
# - _calc_index_counts!(::AbstractExponentsUnbounded, degree::Integer) that makes sure that the cache matrix contains valid
#   entries for at least degree+1 rows
Base.IteratorSize(::Type{<:AbstractExponentsUnbounded}) = Base.IsInfinite()

Base.iterate(e::AbstractExponentsUnbounded{<:Any,I}) where {I<:Integer} =
    ExponentIndices(e, one(I), 0), (one(I), 0, one(I))

function Base.iterate(e::AbstractExponentsUnbounded{<:Any,I}, (index, degree, remainingdeg)::Tuple{I,Int,I}) where {I<:Integer}
    if isone(remainingdeg)
        counts, success = index_counts(e, degree +1)
        @assert(success)
        @inbounds return ExponentIndices(e, index + one(I), degree +1),
            (index + one(I), degree +1, counts[degree+2, 1] - counts[degree+1, 1])
    else
        return ExponentIndices(e, index + one(I), degree), (index + one(I), degree, remainingdeg - one(I))
    end
end

Base.IteratorSize(::Type{<:ExponentsVectorIterator{<:Any,<:AbstractExponentsUnbounded}}) = Base.IsInfinite()

"""
    index_counts(::AbstractExponents{N,I}, degree::Integer) -> Tuple{Matrix{I},Bool}

Safe version of the above. If the boolean in the result is `true`, the matrix will have at least `degree+1` rows, i.e., all
entries up to `degree` are present. If it is false, the requested degree is not present in the exponents and the matrix will
have fewer rows. Note that `true` does not mean that the degree is actually present in the exponents, only that its information
has been calculated.
"""
@inline function index_counts(e::AbstractExponentsUnbounded, degree::Integer)
    # index_counts will be interpreted in a two-dimensional way. The maximal number of variables is fixed, so this will
    # be the fast dimension; the degrees may increase, which leads to appending elements. This will be the slow dimension.
    @boundscheck 0 ≤ degree || throw(ArgumentError("The degree must be nonnegative"))
    @inline(_has_index_counts(e, degree)) || @noinline(_calc_index_counts!(e, degree)) # zero degree = row 1
    # ^ would be good if we could insert an unlikely hint here
    counts = index_counts(unsafe, e)
    return counts, true
end

include("./ExponentsAll.jl")
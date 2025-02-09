export exponents_from_index

struct ExponentIndices{I<:Integer,E<:AbstractExponents{<:Any,I}} <: AbstractVector{Int}
    parent::E
    index::I
    degree::Int
end

"""
    exponents_from_index(unsafe, ::AbstractExponents{N,I}, index::I[, degree::Int])

Unsafe variant of [`exponents_from_index`](@ref exponents_from_index):
assumes that the cache for degree `degree` has already been populated, and that the exponent set contains `index`. If `degree`
is omitted, it is calculated using the unsafe variant of
[`degree_from_index`](@ref degree_from_index(::Unsafe, ::AbstractExponents{N,I}, ::I) where {N,I<:Integer}).
"""
exponents_from_index(::Unsafe, e::AbstractExponents{<:Any,I}, index::I, degree::Int) where {I<:Integer} =
    ExponentIndices(e, index, degree)

# do it manually instead of as a default parameter so that we don't get a missing docstring warning
exponents_from_index(::Unsafe, e::AbstractExponents{<:Any,I}, index::I) where {I<:Integer} =
    exponents_from_index(unsafe, e, index, degree_from_index(unsafe, e, index))

"""
    exponents_from_index(::AbstractExponents{N,I}, index::I[, degree::Int])

Calculates the exponents that are associated with the monomial index `index` in `N` variables within a given exponent set.
The return value will be a lazy implementation of `AbstractVector{Int}` (though iterating is more efficient than indexing).
`degree` must match the degree of the index, if it is specified.

See also [`degree_from_index`](@ref).
"""
exponents_from_index(e::AbstractExponents{<:Any,I}, index::I) where {I<:Integer} =
    exponents_from_index(e, index, degree_from_index(e, index))

Base.IndexStyle(::Type{<:ExponentIndices}) = IndexLinear()
Base.size(::(ExponentIndices{I,<:AbstractExponents{N,I}} where {I<:Integer})) where {N} = (N,)

Base.@assume_effects :foldable :nothrow @inline function Base.getindex(efi::(ExponentIndices{I,<:AbstractExponents{N,I}} where {I<:Integer}),
    varidx::Integer) where {N}
    @boundscheck checkbounds(efi, varidx)
    iter = iterate(efi)::Tuple
    for _ in 2:varidx
        iter = iterate(efi, iter[2])::Tuple
    end
    return iter[1]
end

# Standardized part of the iterator state: degree, index of next variable must be first. Do not change.

# we want the iterator-based versions for copying
Base.@propagate_inbounds Base.copyto!(dest::AbstractArray, src::ExponentIndices) =
    @invoke copyto!(dest::AbstractArray, src::Any)
Base.@propagate_inbounds Base.copyto!(dest::AbstractArray, dstart::Integer, src::ExponentIndices) =
    @invoke copyto!(dest::AbstractArray, dstart::Integer, src::Any)
Base.@propagate_inbounds Base.copyto!(dest::AbstractArray, dstart::Integer, src::ExponentIndices, sstart::Integer) =
    @invoke copyto!(dest::AbstractArray, dstart::Integer, src::Any, sstart::Integer)
Base.@propagate_inbounds Base.copyto!(dest::AbstractArray, dstart::Integer, src::ExponentIndices, sstart::Integer, n::Integer) =
    @invoke copyto!(dest::AbstractArray, dstart::Integer, src::Any, sstart::Integer, n::Integer)

# potentially more efficient shortcut
@inline function exponents_to_index(e::AbstractExponents, exponents::ExponentIndices, ::Int, ::Val{false}=Val(false))
    e === exponents.parent && return exponents.index
    index_counts(e, exponents.degree)
    return convert_index(unsafe, e, exponents.parent, exponents.index)
end

Base.sum(exponents::ExponentIndices; init=0) = exponents.degree
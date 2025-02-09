module FastVector

export FastVec, prepare_push!, unsafe_push!, unsafe_append!, unsafe_prepend!, finish!
using ..PolynomialOptimization: @assert, resizable_array

# Julia exploits the sizehint!/push! combination very poorly (https://github.com/JuliaLang/julia/issues/24909). Here, we define
# a (limited) fast vector, which rectifies this issue. It is similar to the PushVectors.jl package, but allows for even more
# speed by allowing to drop the bounds check on push!. It also uses Julia's precise algorithm to determine the next size when
# growing.
"""
    FastVec{V}(undef, n::Integer; buffer::Integer=n)

Creates a new FastVec, which is a vector of elements of type `V`. The elements are initially undefined; there are `n` items.
The vector has a capacity of size `buffer`; while it can hold arbitrarily many, pushing to the vector is fast as long as the
capacity is not exceeded.
"""
mutable struct FastVec{V} <: AbstractVector{V}
    # We could in principle make this more compact by redefining Julia's internal jl_array_t structure (which already has the
    # length, but let's keep up compatibility)
    const data::Vector{V}
    len::Int

    function FastVec{V}(::UndefInitializer, n::Integer; buffer::Integer=n) where {V}
        buffer ≥ n || error("The buffer must not be smaller than the number of elements")
        new{V}(isbitstype(V) ? resizable_array(V, buffer) : Vector{V}(undef, buffer), n)
    end
end

"""
    FastVec{V}(; buffer::Integer=0)

Shorthand to create an empty FastVec with a certain buffer size.
"""
FastVec{V}(; buffer::Integer=0) where {V} = FastVec{V}(undef, 0; buffer)

Base.size(v::FastVec) = (v.len,)
Base.strides(v::FastVec) = strides(v.data)
Base.elsize(v::FastVec) = Base.elsize(v.data)

iterate(v::FastVec, i=1) = (@inline; (i % UInt) - 1 < length(v) ? (@inbounds v.data[i], i + 1) : nothing)

# Similar to unsafe_cast in SimplePolynomials, but we only want to skip the checks if our given type is not larger than the
# target, because only then we can ensure validity base on the nonnegativity of the length.
unsafe_upcast(T::Type{<:Signed}, x::Signed) = T(x)
unsafe_upcast(T::Type{<:Unsigned}, x::Unsigned) = T(x)
function unsafe_upcast(T::Type{S}, x::Unsigned) where {S<:Signed} # assume the unsigned value does not have the top bit set
    if sizeof(T) == sizeof(x)
        return Core.bitcast(T, x)
    else
        return T(x)
    end
end
function unsafe_upcast(T::Type{U}, x::Signed) where {U<:Unsigned} # assume the signed value is not negative
    if sizeof(T) == sizeof(x)
        return Core.bitcast(T, x)
    else
        return T(x)
    end
end

Base.@propagate_inbounds function Base.getindex(v::FastVec, i::Int)
    @boundscheck checkbounds(v, i)
    @inbounds return v.data[i]
end

Base.IndexStyle(::FastVec) = IndexLinear()

Base.@inline function Base.setindex!(v::FastVec{V}, el, i::Int) where {V}
    @boundscheck checkbounds(v, i)
    @inbounds return v.data[i] = convert(V, el)
end

Base.length(v::FastVec) = v.len

"""
    sizehint!(v::FastVec, len::Integer)

Changes the size of the internal buffer that is kept available to quickly manage pushing into the vector. If len is smaller
than the actual length of this vector, this is a no-op.
"""
@inline function Base.sizehint!(v::FastVec, len::Integer)
    len ≥ length(v.data) && resize!(v.data, unsafe_upcast(UInt, len))
    return v
end

"""
    empty!(v::FastVec)

Clears the vector without freeing the internal buffer.
"""
@inline function Base.empty!(v::FastVec)
    v.len = 0
    return v
end

"""
    prepare_push!(v::FastVec, new_items::Integer)

Prepares pushing (or appending) at last `new_items` in the future, in one or multiple calls. This ensures that the internal
buffer is large enough to hold all the new items that will be pushed without allocations in between.
"""
@inline function prepare_push!(v::FastVec, new_items::Integer)
    @boundscheck(new_items ≥ 0 || error("prepare_push! expects a positive number of items"))
    length(v.data) - v.len < new_items && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len + new_items)))
    return v
end

"""
    push!(v::FastVec, el)

Adds `el` to the end of `v`, increasing the length of `v` by one. If there is not enough space, grows `v`.
Note that if you made sure beforehand that the capacity of `v` is sufficient for the addition of this element, consider calling
[`unsafe_push!`](@ref) instead, which avoids the length check.
"""
@inline function Base.push!(v::FastVec{V}, el) where {V}
    elV = convert(V, el)
    v.len += 1
    length(v.data) < v.len && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len)))
    @inbounds v.data[v.len] = elV
    return v
end

"""
    unsafe_push!(v::FastVec, el)

Adds `el` to the end of `v`, increasing the length of `v` by one.
This function assumes that the internal buffer of `v` holds enough space to add at least one element; if this is not the case,
it will lead to memory corruption. Call [`push!`](@ref) instead if you cannot guarantee the necessary buffer size.
"""
@inline function unsafe_push!(v::FastVec{V}, el) where {V}
    elV = convert(V, el)
    v.len += 1
    @assert(length(v.data) ≥ v.len)
    @inbounds v.data[v.len] = elV
    return v
end

@inline function unsafe_push!(v::FastVec{V}, els...) where {V}
    # we don't need calls to copyto! here - the size of els is completely known and probably quite small, so let's inline
    @assert(length(v.data) ≥ v.len + length(els))
    for el in els
        @inbounds v.data[v.len += 1] = convert(V, el)
    end
    return v
end

"""
    insert!(v::FastVec, index::Integer, el)

Insert an `el` into `v` at the given `index`. `index` is the index of `item` in the resulting `v`.
Note that if you made sure beforehand that the capacity of `v` is sufficient for the addition of this element, consider calling
[`unsafe_insert!`](@ref) instead, which avoids the length check.
"""
@inline function Base.insert!(v::FastVec{V}, index::Integer, el) where {V}
    @boundscheck(1 ≤ index ≤ v.len +1 || throw(BoundsError(v, index)))
    elV = convert(V, el)
    v.len += 1
    length(v.data) < v.len && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len)))
    @inbounds copyto!(v.data, index +1, v.data, index, v.len - index)
    @inbounds v.data[index] = el
    return v
end

"""
    unsafe_insert!(v::FastVec, index::Integer, el)

Insert an `el` into `v` at the given `index`. `index` is the index of `item` in the resulting `v`.
This function assumes that the internal buffer of `v` holds enough space to add at least one element; if this is not the case,
it will lead to memory corruption. Call [`insert!`](@ref) instead if you cannot guarantee the necessary buffer size.
"""
@inline function unsafe_insert!(v::FastVec{V}, index::Integer, el) where {V}
    @assert(1 ≤ index ≤ v.len +1 || throw(BoundsError(v, index)))
    elV = convert(V, el)
    v.len += 1
    @inbounds copyto!(v.data, index +1, v.data, index, v.len - index)
    @inbounds v.data[index] = el
    return v
end

"""
    append!(v::FastVec, els)

Appends all items in `els` to the end of `v`, increasing the length of `v` by `length(els)`. If there is not enough space,
grows `v`.
Note that if you made sure beforehand that the capacity of `v` is sufficient for the addition of these elements, consider
calling [`unsafe_append!`](@ref) instead, which avoids the length check.
"""
@inline function Base.append!(v::FastVec, els)
    oldlen = v.len
    v.len += length(els)
    length(v.data) < v.len && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len)))
    @inbounds copyto!(v.data, oldlen +1, els, 1, length(els))
    return v
end

@inline function Base.append!(v::FastVec{V}, els::FastVec{V}) where {V}
    oldlen = v.len
    v.len += els.len
    length(v.data) < v.len && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len)))
    @inbounds copyto!(v.data, oldlen +1, els.data, 1, els.len)
    return v
end

"""
    unsafe_append!(v::FastVec, els)

Appends all items in `els` to the end of `v`, increasing the length of `v` by `length(els)`.
This function assumes that the internal buffer of `v` holds enough space to add at least all elements in `els`; if this is not
the case, it will lead to memory corruption. Call [`append!`](@ref) instead if you cannot guarantee the necessary buffer size.
"""
@inline function unsafe_append!(v::FastVec, els)
    oldlen = v.len
    v.len += length(els)
    @assert(length(v.data) ≥ v.len)
    @inbounds copyto!(v.data, oldlen +1, els, 1, length(els))
    return v
end

@inline function unsafe_append!(v::FastVec{V}, els::FastVec{V}) where {V}
    oldlen = v.len
    v.len += els.len
    @assert(length(v.data) ≥ v.len)
    @inbounds copyto!(v.data, oldlen +1, els.data, 1, els.len)
    return v
end

"""
    prepend!(v::FastVec, els)

Prepends all items in `els` to the beginning of `v`, increasing the length `v` by `length(els)`. If there is not enough space,
grows `v`.
Note that if you made sure beforehand that the capacity of `v` is sufficient for the addition of these elements, consider
calling [`unsafe_prepend!`](@ref) instead, which avoids the length check.
"""
@inline function Base.prepend!(v::FastVec, els)
    oldlen = v.len
    v.len += length(els)
    length(v.data) < v.len && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len)))
    @inbounds copyto!(v.data, length(els) +1, v.data, 1, oldlen) # must always work correctly with overlapping
    @inbounds copyto!(v.data, els)
    return v
end

@inline function Base.prepend!(v::FastVec{V}, els::FastVec{V}) where {V}
    oldlen = v.len
    v.len += els.len
    length(v.data) < v.len && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len)))
    @inbounds copyto!(v.data, els.len +1, v.data, 1, oldlen) # must always work correctly with overlapping
    @inbounds copyto!(v.data, 1, els.data, 1, els.len)
    return v
end

"""
    unsafe_prepend!(v::FastVec, els::AbstractVector)

Prepends all items in `els` to the beginning of `v`, increasing the length `v` by `length(els)`.
This function assumes that the internal buffer of `v` holds enough space to add at least all elements in `els`; if this is not
the case, it will lead to memory corruption. Call [`prepend!`](@ref) instead if you cannot guarantee the necessary buffer size.
"""
@inline function unsafe_prepend!(v::FastVec, els)
    oldlen = v.len
    v.len += length(els)
    @assert(length(v.data) ≥ v.len)
    @inbounds copyto!(v.data, length(els) +1, v.data, 1, oldlen) # must always work correctly with overlapping
    @inbounds copyto!(v.data, els)
    return v
end

@inline function unsafe_prepend!(v::FastVec{V}, els::FastVec{V}) where {V}
    oldlen = v.len
    v.len += els.len
    @assert(length(v.data) ≥ v.len)
    @inbounds copyto!(v.data, els.len +1, v.data, 1, oldlen) # must always work correctly with overlapping
    @inbounds copyto!(v.data, 1, els.data, 1, els.len)
    return v
end

"""
    splice!(v::FastVec, index::Integer, [replacement]) -> item

Remove the item at the given index, and return the removed item. Subsequent items are shifted left to fill the resulting gap.
If specified, replacement values from an ordered collection will be spliced in place of the removed item.
No unsafe version of this function exists.
To insert `replacement` before an index `n` without removing any items, use `splice!(collection, n:n-1, replacement)`.
"""
function Base.splice!(v::FastVec, i::Integer, ins=Base._default_splice)
    x = v[i]
    m = length(ins)
    @inbounds if m == 0
        deleteat!(v, i)
    elseif m == 1
        v.data[i] = first(ins)
    else
        v.len += m -1
        length(v.data) < v.len && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len)))
        copyto!(v.data, i + m, v.data, i +1, v.len - m - i +1)
        k = 0
        for x in ins
            v.data[i+k] = x
            k += 1
        end
    end
    return x
end

function Base.splice!(v::FastVec, r::AbstractUnitRange{<:Integer}, ins=Base._default_splice)
    x = v[r]
    m = length(ins)
    if m == 0
        deleteat!(v, r)
        return x
    end

    n = v.len
    f = first(r)
    l = last(r)
    d = length(r)

    @inbounds if m < d
        delta = d - m
        copyto!(v.data, l - delta +1, v.data, l +1, n - f - delta)
        v.len -= delta
    elseif m > d
        delta = m - d
        v.len += delta
        length(v.data) < v.len && resize!(v.data, overallocation(unsafe_upcast(UInt, v.len)))
        copyto!(v.data, l + delta +1, v.data, l +1, n - l)
    end

    k = f
    @inbounds for y in ins
        v.data[k] = y
        k += 1
    end
    return x
end

"""
    similar(v::FastVec)

Creates a FastVec of the same type and length as `v`. All the common variants to supply different element types or lengths are
also available; when changing the length, you might additionally specify the keyword argument `buffer` that also allows to
change the internal buffer size.
"""
Base.similar(v::FastVec{V}) where {V} = FastVec{V}(undef, v.len, buffer=length(v.data))
Base.similar(v::FastVec{V}, ::Type{S}) where {V,S} = FastVec{S}(undef, v.len, buffer=length(v.data))
Base.similar(::FastVec{V}, len::Integer; buffer::Integer=overallocation(len)) where {V} = FastVec{V}(undef, len; buffer)
Base.similar(::FastVec{V}, ::Type{S}, len::Integer; buffer::Integer=overallocation(len)) where {V,S} =
    FastVec{S}(undef, len; buffer)

Base.unsafe_convert(::Type{Ptr{V}}, v::FastVec{V}) where {V} = Base.unsafe_convert(Ptr{V}, v.data)

"""
    copyto!(dest::FastVec, doffs::Integer, src::FastVec, soffs::Integer, n::Integer)
    copyto!(dest::Array, doffs::Integer, src::FastVec, soffs::Integer, n::Integer)
    copyto!(dest::FastVec, doffs::Integer, src::Array, soffs::Integer, n::Integer)

Implements the standard `copyto!` operation between `FastVec`s and also mixed with source or destination as an array.
"""
@inline function Base.copyto!(dest::FastVec, doffs::Integer, src::FastVec, soffs::Integer, n::Integer)
    n == 0 && return dest
    n > 0 || Base._throw_argerror()
    @boundscheck checkbounds(dest, doffs:doffs+n-1)
    @boundscheck checkbounds(src, soffs:soffs+n-1)
    unsafe_copyto!(dest.data, doffs, src.data, soffs, n)
    return dest
end

@inline function Base.copyto!(dest::Array, doffs::Integer, src::FastVec, soffs::Integer, n::Integer)
    n == 0 && return dest
    n > 0 || Base._throw_argerror()
    @boundscheck checkbounds(dest, doffs:doffs+n-1)
    @boundscheck checkbounds(src, soffs:soffs+n-1)
    unsafe_copyto!(dest, doffs, src.data, soffs, n)
    return dest
end

@inline function Base.copyto!(dest::FastVec, doffs::Integer, src::Array, soffs::Integer, n::Integer)
    n == 0 && return dest
    n > 0 || Base._throw_argerror()
    @boundscheck checkbounds(dest, doffs:doffs+n-1)
    @boundscheck checkbounds(src, soffs:soffs+n-1)
    unsafe_copyto!(dest.data, doffs, src, soffs, n)
    return dest
end

"""
    resize!(v::FastVec, n::Integer)

Ensures that the internal buffer can hold at least `n` items (meaning that larger buffers will not be shrunk, but smaller ones
will be increased to exactly `n`) and sets the length of the vector to `n`.
"""
function Base.resize!(v::FastVec, n::Integer)
    n > length(v.data) && resize!(v.data, unsafe_upcast(UInt, n)) # assume we know what we are doing, so no overallocation here
    v.len = n
    return v
end

"""
    deleteat!(a::FastVec, i)

Remove the item at the given `i` and return the modified `a`. Subsequent items are shifted to fill the resulting gap. The
internal buffer size is not modified. The index must be either an integer or a unit range.
"""
function Base.deleteat!(v::FastVec, i::Integer)
    @inbounds for j in i:v.len-1
        v.data[j] = v.data[j+1]
    end
    v.len -= 1
    return v
end

Base.@propagate_inbounds function Base.deleteat!(v::FastVec, i::AbstractUnitRange)
    copyto!(v.data, first(i), v.data, last(i) +1, v.len - last(i))
    v.len -= length(i)
    return v
end

@inline function Base._deleteend!(v::FastVec, delta::Integer)
    @boundscheck checkbounds(v, v.len - delta +1)
    v.len -= delta
    return v
end

function Base._growbeg!(v::FastVec, delta::Integer)
    delta = Int(delta)
    delta == 0 && return # avoid attempting to index off the end
    delta >= 0 || throw(ArgumentError("grow requires delta >= 0"))
    if v.len + delta ≤ length(v.data)
        @inbounds copyto!(v.data, delta +1, v.data, 1, v.len)
    else
        newdelta = (v.len + delta) - length(v.data)
        Base._growbeg!(v.data, newdelta)
        @inbounds copyto!(v.data, delta +1, v.data, newdelta +1, v.len) # unfortunately, we have to copy twice
    end
    v.len += delta
    return v
end

function Base._growend!(v::FastVec, delta::Integer)
    delta = Int(delta)
    delta >= 0 || throw(ArgumentError("grow requires delta >= 0"))
    v.len + delta ≤ length(v.data) || Base._growend!(v.data, delta - (length(v.data) - v.len))
    v.len += delta
    return v
end

function Base._growat!(v::FastVec, i::Integer, delta::Integer)
    delta = Int(delta)
    i = Int(i)
    i == 1 && return Base._growbeg!(v, delta)
    len = v.len
    i == len + 1 && return Base._growend!(v, delta)
    delta ≥ 0 || throw(ArgumentError("grow requires delta >= 0"))
    1 < i ≤ len || throw(BoundsError(v, i))
    if v.len + delta ≤ length(v.data)
        @inbounds copyto!(v.data, i + delta, v.data, i, v.len - i +1)
    else
        newdelta = (v.len + delta) - length(v.data)
        Base._growat!(v.data, i, newdelta)
        @inbounds copyto!(v.data, i + delta, v.data, i + newdelta, v.len - i +1)
    end
    v.len += delta
    return v
end

# inspired by KristofferC's PushVector
"""
    finish!(v::FastVec)

Returns the `Vector` representation that internally underlies the FastVec instance. This function should only be called when no
further operations on the FastVec itself are carried out.
"""
finish!(v::FastVec) = resize!(v.data, unsafe_upcast(UInt, v.len))

# This is how Julia internally grows vectors:
# size_t overallocation(size_t maxsize)
# {
#     if (maxsize < 8)
#         return 8;
#     // compute maxsize = maxsize + 4*maxsize^(7/8) + maxsize/8
#     // for small n, we grow faster than O(n)
#     // for large n, we grow at O(n/8)
#     // and as we reach O(memory) for memory>>1MB,
#     // this means we end by adding about 10% of memory each time
#     int exp2 = sizeof(maxsize) * 8 -
# #ifdef _P64
#         __builtin_clzll(maxsize);
# #else
#         __builtin_clz(maxsize);
# #endif
#     maxsize += ((size_t)1 << (exp2 * 7 / 8)) * 4 + maxsize / 8;
#     return maxsize;
# }
function overallocation(maxsize::I) where {I<:Integer}
    maxsize < I(8) && return I(8)
    exp2 = I(8) * sizeof(maxsize) - leading_zeros(maxsize)
    return maxsize + I(4(1 << (7exp2 ÷ 8)) + maxsize ÷ 8)
end

end
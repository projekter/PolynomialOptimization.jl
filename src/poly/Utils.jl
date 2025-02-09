export ConstantVector

# Define a Missing replacement, as we want to provide some conveniences that probably would not be welcome in Base...
const XorTX{X} = Union{X,<:Type{X}}

# Find the smallest possible type that can hold values between 0 and maxval
Base.@assume_effects :total function smallest_unsigned(maxval::Integer)
    maxval ≤ typemax(UInt8) && return UInt8
    maxval ≤ typemax(UInt16) && return UInt16
    maxval ≤ typemax(UInt32) && return UInt32
    return UInt64
end
smallest_unsigned(::Val{T}) where {T} = Val(typeintersect(T, Unsigned))
smallest_unsigned(::T) where {T} = typeintersect(T, Unsigned)

function _sortedallunique(v::AbstractVector)
    for (x, y) in zip(v, Iterators.drop(v, 1))
        x == y && return false
    end
    return true
end

struct ConstantVector{T} <: AbstractVector{T}
    value::T
    length::Int
end

Base.IndexStyle(::Type{<:ConstantVector}) = IndexLinear()
Base.size(cv::ConstantVector) = (cv.length,)
@inline function Base.getindex(cv::ConstantVector, i)
    @boundscheck checkbounds(cv, i)
    return cv.value
end
Base.iterate(cv::ConstantVector, rem=cv.length) = rem ≤ 0 ? nothing : (cv.value, rem -1)
@inline function elementwise(f, a::ConstantVector, b::ConstantVector)
    @assert(length(a) == length(b))
    return ConstantVector(f(a.value, b.value), a.length)
end
@inline function elementwise(f, a::ConstantVector, b)
    @assert(length(a) == length(b))
    return f.((a.value,), b)
end
@inline function elementwise(f, a, b::ConstantVector)
    @assert(length(a) == length(b))
    return f.(a, (b.value,))
end
Base.@propagate_inbounds elementwise(f, a, b) = f.(a, b)
function elementwise!(::ConstantVector, f, a, b)
    @assert(length(into) == length(a))
    return elementwise(f, a, b)
end
Base.@propagate_inbounds elementwise!(into::AbstractVector, f, a::ConstantVector, b) = into .= f.((a.value,), b)
Base.@propagate_inbounds elementwise!(into::AbstractVector, f, a, b::ConstantVector) = into .= f.(a, (b.value,))
Base.@propagate_inbounds elementwise!(into::AbstractVector, f, a, b) = into .= f.(a, b)

function intersect_sorted(a, b)
    result = FastVec{promote_type(eltype(a), eltype(b))}(buffer=min(length(a), length(b)))
    ia = iterate(a)
    ib = iterate(b)
    (isnothing(ia) || isnothing(ib)) && return finish!(result)
    while true
        if ia[1] == ib[1]
            unsafe_push!(result, ia[1])
            (isnothing((ia = iterate(a, ia[2]);)) ||
                isnothing((ib = iterate(b, ib[2]);))) && return finish!(result)
        elseif ia[1] > ib[1]
            isnothing((ib = iterate(b, ib[2]);)) && return finish!(result)
        else
            isnothing((ia = iterate(a, ia[2]);)) && return finish!(result)
        end
    end
end
intersect_sorted(a::AbstractVector, b::AbstractUnitRange) =
    @view(a[searchsortedfirst(a, first(b)):searchsortedlast(a, last(b))])
intersect_sorted(a::AbstractUnitRange, b::AbstractVector) = intersect_sorted(b, a)
intersect_sorted(a::AbstractUnitRange, b::AbstractUnitRange) = intersect(a, b)
# TODO (maybe): intersect_sorted for two vectors, using binary or exponential search skipping...

unsafe_cast(T::Type{<:Signed}, x::Signed) = T(x)
unsafe_cast(T::Type{<:Unsigned}, x::Unsigned) = T(x)
function unsafe_cast(T::Type{S}, x::Unsigned) where {S<:Signed} # assume the unsigned value does not have the top bit set
    if sizeof(T) == sizeof(x)
        return Core.bitcast(T, x)
    elseif sizeof(T) < sizeof(x)
        return Core.trunc_int(T, x)
    else
        return T(x)
    end
end
function unsafe_cast(T::Type{U}, x::Signed) where {U<:Unsigned} # assume the signed value is not negative
    if sizeof(T) == sizeof(x)
        return Core.bitcast(T, x)
    elseif sizeof(T) < sizeof(x)
        return Core.trunc_int(T, x)
    else
        return T(x)
    end
end
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

mutable struct ConstantVector{T} <: AbstractVector{T}
    length::Int
    value::T

    ConstantVector{T}(value::T, length::Integer) where {T} = new{T}(length, value)
    ConstantVector{T}(::UndefInitializer, length::Integer) where {T} = new{T}(length)
end

ConstantVector(value::T, length::Integer) where {T} = ConstantVector{T}(value, length)

Base.IndexStyle(::Type{<:ConstantVector}) = IndexLinear()
Base.size(cv::ConstantVector) = (cv.length,)
@inline function Base.getindex(cv::ConstantVector, i::Integer)
    @boundscheck checkbounds(cv, i)
    return cv.value
end
@inline function Base.getindex(cv::ConstantVector{T}, i::Union{<:AbstractRange,Colon,<:AbstractArray}) where {T}
    @boundscheck checkbounds(cv, i)
    return ConstantVector{T}(cv.value, length(i))
end
Base.iterate(cv::ConstantVector, rem=cv.length) = rem ≤ 0 ? nothing : (cv.value, rem -1)
function Base._mapreduce(f, op, ::Base.IndexCartesian, a::ConstantVector)
    T = eltype(a)
    isempty(a) && return Base.mapreduce_empty(f, op, T)
    return _mapreducerep(f, op, a.value, a.length)
end
Base._any(f, a::ConstantVector, ::Colon) = !isempty(a) && f(a.value)
Base._all(f, a::ConstantVector, ::Colon) = isempty(a) || f(a.value)
Base.map(f, a::ConstantVector) = ConstantVector(f(a.value), a.length)
_mapreducerep(f, op::Union{typeof(Base.add_sum),typeof(+)}, v, num::Integer) = f(v) * num
_mapreducerep(f, op::Union{typeof(Base.mul_prod),typeof(*)}, v, num::Integer) = f(v) ^ num
_mapreducerep(f, op::Union{typeof(min),typeof(max)}, v, ::Integer) = f(v)
_mapreducerep(f::Base.ExtremaMap, op::typeof(Base._extrema_rf), v, ::Integer) = f(v)
function _mapreducerep(f, op, v, num::Integer)
    val = f(v)
    result = val
    for i in 2:num
        newval = op(result, v)
        isequal(newval, result) && return result # fixed point
        result = newval
    end
    return result
end
Base.sort(a::ConstantVector; kwargs...) = a
Base.sort!(a::ConstantVector; kwargs...) = a
Base.BroadcastStyle(::Type{<:ConstantVector}) = Broadcast.ArrayStyle{ConstantVector}()
Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{ConstantVector}}, ::Type{T}) where {T} =
    ConstantVector{T}(undef, length(bc))
@inline function Base.copyto!(dest::ConstantVector, bc::Broadcast.Broadcasted{Nothing})
    axes(dest) == axes(bc) || Broadcast.throwdm(axes(dest), axes(bc))
    if bc.f === identity && bc.args isa Tuple{ConstantVector}
        A = bc.args[1]
        if axes(dest) == axes(A)
            dest.value = A.value
            return dest
        end
    end
    bc′ = Broadcast.preprocess(dest, bc)
    dest.value = first(bc′)
    return dest
end
# we have to undo the "performance" optimization
@inline function Base.copyto!(dest::ConstantVector, bc::Broadcast.Broadcasted{<:Broadcast.AbstractArrayStyle{0}})
    if bc.f === identity && bc.args isa Tuple{Any} && Broadcast.isflat(bc)
        dest.value = bc.args[1][]
        return dest
    else
        return copyto!(dest, convert(Broadcast.Broadcasted{Nothing}, bc))
    end
end
# ConstantVector style loses to any but scalars
Base.BroadcastStyle(::Broadcast.ArrayStyle{ConstantVector}, ::Broadcast.AbstractArrayStyle{0}) =
    Broadcast.ArrayStyle{ConstantVector}()
Base.BroadcastStyle(::Broadcast.ArrayStyle{ConstantVector}, ::Broadcast.DefaultArrayStyle{0}) =
    Broadcast.ArrayStyle{ConstantVector}()
Base.BroadcastStyle(::Broadcast.ArrayStyle{ConstantVector}, a::Broadcast.DefaultArrayStyle) = a

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
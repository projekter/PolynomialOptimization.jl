# We don't want to import StaticArrays for a simple stack-allocated tuple whose address we can get
mutable struct StackVec{N,T} <: AbstractVector{T}
    data::NTuple{N,T}

    StackVec(x::T...) where {T} = new{length(x),T}(x)

    StackVec{T}() where {T} = new{0,T}(())
end

Base.size(::StackVec{N}) where {N} = (N,)
Base.length(::StackVec{N}) where {N} = N
Base.eltype(::StackVec{<:Any,T}) where {T} = T

@inline function Base.unsafe_convert(::Type{Ptr{T}}, v::StackVec{<:Any,T}) where {T}
    isbitstype(T) || error("StackVec pointer conversion requires bitstype elements")
    return Base.unsafe_convert(Ptr{T}, pointer_from_objref(v))
end

Base.convert(::Type{Vector}, v::StackVec) = collect(v.data)
Base.convert(::Type{Vector{T}}, v::StackVec{<:Any,T}) where {T} = collect(v.data)

@inline function Base.getindex(v::StackVec, i::Integer)
    @boundscheck checkbounds(v, i)
    @inbounds return v.data[i]
end

# not used anywhere - we could theoretically make data constant
#=@inline function Base.setindex!(v::StackVec{N,T} where {N}, i::Integer, x::T) where {T}
    isbitstype(T) || error("StackVec mutation requires bitstype elements")
    @boundscheck checkbounds(v.data, i)
    GC.@preserve v unsafe_store!(Base.unsafe_convert(Ptr{T}, pointer_from_objref(v)), convert(T, x), i)
    return x
end=#
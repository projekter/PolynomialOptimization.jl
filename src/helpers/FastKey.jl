struct FastKey{T<:Integer}
    value::T
end

Base.:(==)(x::FastKey, y::FastKey) = x.value == y.value
Base.hash(x::FastKey, h::UInt) = UInt(x.value) ⊻ h
Base.convert(::Type{T}, x::FastKey) where {T<:Integer} = convert(T, x.value)::T
Base.isless(a::FastKey, b::FastKey) = a.value < b.value
Base.isless(a::FastKey, b) = a.value < b
Base.isless(a, b::FastKey) = a < b.value
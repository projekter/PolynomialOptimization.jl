#region Fixed-length byte string
struct _TupleVector{N,T} <: AbstractVector{T}
    data::NTuple{N,T}
end

Base.IndexStyle(::Type{<:_TupleVector}) = Base.IndexLinear()
Base.size(::_TupleVector{N}) where {N} = (N,)
Base.@propagate_inbounds Base.getindex(t::_TupleVector, i) = t.data[i]

struct FortranString{Len} <: AbstractString
    data::_TupleVector{Len,UInt8}

    function FortranString{Len}(str::AbstractString) where {Len}
        s = sizeof(codeunit(str))
        strlen = s * ncodeunits(str)
        # Fill our data with (little endian) bytes that make up the string. We don't do any conversion, and our output will not
        # be interpreted in any kind of encoding.
        if strlen ≥ Len
            new{Int(Len)}(_TupleVector(ntuple(let s=s
                i -> UInt8((codeunit(str, div(i, s, RoundUp)) >> 8((i -1) % s)) & 0xFF)
            end, Val(Len))))
        else
            new{Int(Len)}(_TupleVector(ntuple(let s=s, strlen=strlen
                i -> i ≤ strlen ? UInt8((codeunit(str, div(i, s, RoundUp)) >> 8((i -1) % s)) & 0xFF) :
                                  (iszero((i - strlen -1) % s) ? 0x20 : 0x00)
            end, Val(Len))))
        end
    end
end

Base.codeunit(::FortranString) = UInt8
Base.codeunit(f::FortranString, i::Integer) = f.data[i]
Base.ncodeunits(::FortranString{Len}) where {Len} = Len
Base.isvalid(f::FortranString, i::Int) = checkbounds(Bool, f, i)
function Base.iterate(f::FortranString, args...)
    iter = iterate(f.data, args...)
    isnothing(iter) && return nothing
    return FortranChar(iter[1]), iter[2]
end
Base.@propagate_inbounds Base.getindex(f::FortranString, i) = FortranChar(f[i])

struct FortranChar <: AbstractChar
    c::UInt8
end

Base.codepoint(c::FortranChar) = UInt32(c.c) # no encoding, we don't care
#endregion

#region 4-byte boolean
struct FortranBool <: Integer
    value::Cint

    FortranBool(value::Bool) = new(Cint(value))
end

Base.convert(::Type{FortranBool}, v::Bool) = FortranBool(v)
Base.convert(::Type{Bool}, v::FortranBool) = !iszero(v.value)
Base.show(io::IO, x::FortranBool) = show(io, !iszero(x.value))
#endregion

#region Fortran Array Descriptors
# Note that LANCELOT is the only part of GALAHAD that does not have a C interface yet. Hence, we first need to provide some
# Fortran magic - unfortunately, all the array arguments are passed as assumed-shape, so that we have to pass array
# descriptors. Warning: Here, we assume that GALAHAD was compiled using GFortran with a version ≥ 9. This is a moderately safe
# assumption, as the Intel compiler does not work out-of-the-box... but in case you wish to adapt this to Intel, this is pretty
# easy, as their array descriptor is well-documented:
# https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-1/handle-fortran-array-descriptors.html
# In contrast, the GFortran array descriptor is not documented properly, the following structure was extracted directly from
# the source code:
# https://github.com/gcc-mirror/gcc/blob/386df7ce7b38ef00e28080a779ef2dfd6949cf15/libgfortran/libgfortran.h#L433-L440
# Note: This changed with GCC 9 (PR37577 and PR34640), when the step to Fortran 2008 was made. Before, there was no dedicated
# DType struct, but a bitmask instead.
struct GFortranArrayDescriptorDType
    elem_len::Csize_t
    version::Cint
    rank::Cchar
    type::Cchar
    attribute::Cshort
end

struct GFortranArrayDescriptorDimension
    stride::Cptrdiff_t
    lower_bound::Cptrdiff_t
    upper_bound::Cptrdiff_t
end

struct GFortranArrayDescriptor{T,N} <: AbstractArray{T,N}
    data::Ptr{T} # cannot be Ref, as this may be C_NULL
    offset::Cptrdiff_t
    dtype::GFortranArrayDescriptorDType
    span::Cptrdiff_t
    dim::NTuple{N,GFortranArrayDescriptorDimension}

    GFortranArrayDescriptor{T,N}(::Missing=missing, _=missing) where {T,N} =
        new{T,N}(Ptr{T}(C_NULL)) # the rest is just uninitialized
    GFortranArrayDescriptor{T,N}(data::Ptr{T}, offset::Cptrdiff_t, dtype::GFortranArrayDescriptorDType, span::Cptrdiff_t,
        dim::NTuple{N,GFortranArrayDescriptorDimension}) where {T,N} = new{T,N}(data, offset, dtype, span, dim)
    GFortranArrayDescriptor{T,N}(a::AbstractArray{T,N}, pushref::AbstractVector) where {T,N} =
        convert(GFortranArrayDescriptor{T,N}, a, pushref)
end

gfortran_array_element_type(::Type{<:Integer}) = Cchar(1)
gfortran_array_element_type(::Type{FortranBool}) = Cchar(2)
gfortran_array_element_type(::Type{<:Real}) = Cchar(3)
gfortran_array_element_type(::Type{<:Complex}) = Cchar(4)
gfortran_array_element_type(::Any) = Cchar(5)
gfortran_array_element_type(::Type{<:Union{<:AbstractChar,AbstractString}}) = Cchar(6)
# BT_CLASS, BT_PROCEDURE, BT_HOLLERITH, BT_VOID, BT_ASSUMED, BT_UNION, BT_BOZ not supported

Base.convert(G::Type{<:GFortranArrayDescriptor}, X::AbstractArray{T,N}, args...) where {T,N} =
    convert(GFortranArrayDescriptor{T,N}, X, args...)
function Base.convert(TT::Type{GFortranArrayDescriptor{T,N}}, X::AbstractArray{T,N}, pushref::AbstractVector) where {T,N}
    N ≤ 15 || throw(MethodError(convert, (TT, X)))
    s = strides(X)
    a = axes(X)
    eltype(a) <: AbstractUnitRange || throw(MethodError(convert, (TT, X)))
    # Careful - GFortran assigns C_NULL to data if it is not allocated. To distinguish this from an allocated zero-length
    # array, a one-byte malloc call is performed for empty arrays. To mimick this, we must at some place store a reference to
    # this malloc data. Note that as we directly alias the Julia data, we won't allow the Fortran code to alter our data
    # anyway, so we can just create a Ref. But we need to store this somewhere for bookkeeping.
    if isempty(X)
        data = Ptr{T}(Base.unsafe_convert(Ptr{UInt8}, push!(pushref, Ref{UInt8}(0))))
    else
        push!(pushref, X)
        data = pointer(X)
    end
    return TT(
        data,
        -sum(first(aᵢ) * sᵢ for (aᵢ, sᵢ) in zip(a, s); init=0),
        GFortranArrayDescriptorDType(
            sizeof(T),
            0, # version appears to be unused
            N,
            gfortran_array_element_type(T),
            0 # attributes seem to be unused - in particular, GFC_ARRAY_ASSUMED_SHAPE_CONT or GFC_ARRAY_ASSUMED_SHAPE appear to
              # be irrelevant
        ),
        sizeof(T), # span appears to be identical to the elem_len (at least in all considered cases)
        ntuple(let s=s, a=a; i -> @inbounds GFortranArrayDescriptorDimension(s[i], first(a[i]), last(a[i])) end, N),
    )
end
Base.convert(::Type{Array}, D::GFortranArrayDescriptor{T,N}) where {T,N} = convert(Array{T,N}, D)
function Base.convert(::Type{Array{T,N}}, D::GFortranArrayDescriptor{T,N}) where {T,N}
    D.data == C_NULL && return Array{T,N}(undef, ntuple(_ -> 0, Val(N)))
    @assert(sizeof(T) == D.dtype.elem_len && N == D.dtype.rank)
    @assert(isone(D.dim[end].stride) &&
        all(i -> D.dim[i].stride == D.dim[i+1].stride * (D.dim[i+1].upper_bound - D.dim[i+1].lower_bound +1), 1:N-1))
    # We only do an imperfect translation, as our arrays will always start at 1
    return unsafe_wrap(Array, D.data, NTuple{N,Int}(d.upper_bound - d.lower_bound + 1 for d in D.dim))
end
const FVector{T} = GFortranArrayDescriptor{T,1}
const FMatrix{T} = GFortranArrayDescriptor{T,2}

Base.size(g::GFortranArrayDescriptor{<:Any,N}) where {N} =
    g.data == C_NULL ? ntuple(_ -> 0, Val(N)) : ntuple(i -> Int(g.dim[i].upper_bound - g.dim[i].lower_bound +1), Val(N))
# For efficiency reasons, we also use sizeof(T) here instead of g.dtype.elem_len, which should really be the same. This is
# still not as efficient as a Julia Vector, as we cannot statically infer anything about contiguity.
@inline @generated function Base.getindex(g::GFortranArrayDescriptor{T,N}, I::Vararg{Int,N}) where {T,N}
    addr = Expr(:call, :+, :(g.offset))
    for i in 1:N
        push!(addr.args, :(g.dim[$i].stride * I[$i]))
    end
    quote
        @boundscheck checkbounds(g, CartesianIndex(I))
        return unsafe_load(g.data + sizeof(T) * $addr)
    end
end
@inline @generated function Base.setindex!(g::GFortranArrayDescriptor{T,N}, v, I::Vararg{Int,N}) where {T,N}
    addr = Expr(:call, :+, :(g.offset))
    for i in 1:N
        push!(addr.args, :(g.dim[$i].stride * I[$i]))
    end
    quote
        @boundscheck checkbounds(g, CartesianIndex(I))
        return unsafe_store!(g.data + sizeof(T) * $addr, v)
    end
end
Base.axes(g::GFortranArrayDescriptor{<:Any,N}) where {N} =
    g.data == C_NULL ? ntuple(_ -> 0, Val(N)) : ntuple(i -> g.dim[i].lower_bound:g.dim[i].upper_bound, Val(N))
Base.strides(g::GFortranArrayDescriptor{<:Any,N}) where {N} =
    g.data == C_NULL ? ntuple(_ -> 0, Val(N)) : ntuple(i -> Int(g.dim[i].stride))
#endregion
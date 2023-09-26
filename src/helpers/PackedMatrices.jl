module PackedMatrices

export PackedMatrix, PackedMatrixUpper, PackedMatrixLower, packed_isupper, packed_islower, packed_isscaled, packedsize,
    rmul_diags!, rmul_offdiags!, lmul_diags!, lmul_offdiags!, packed_scale!, packed_unscale!, spev!, spevd!, spevx!, pptrf!,
    spmv!, spr!, tpttr!, trttp!, gemmt!, eigmin!, eigmax!

using LinearAlgebra, SparseArrays

# Symmetric square matrix. Depending on Fmt, we only store the upper or lower triangle in col-major vectorized or scaled
# vectorized form.
# The PackedMatrix can be effciently broadcast to other matrices, which works on the vector representation. It can also be
# receive values from generic matrices by copyto! or broadcasting. Using it in a mixed chain of broadcastings of different type
# is not implemented and will potentially lead to logical errors (in the sense that the types will not match) or even segfaults
# (as the correct index mapping is not implemented).
struct PackedMatrix{R,V<:AbstractVector{R},Fmt} <: AbstractMatrix{R}
    dim::Int
    data::V

    function PackedMatrix(dim::Integer, data::AbstractVector{R}, type::Symbol=:U) where {R}
        chkpacked(dim, data)
        type ∈ (:U, :L, :US, :LS) || error("Type must be one of :U, :L, :US, or :LS")
        return new{R,typeof(data),type}(dim, data)
    end

    function PackedMatrix{R}(::UndefInitializer, dim::Integer, type::Symbol=:U) where {R}
        type ∈ (:U, :L, :US, :LS) || error("Type must be one of :U, :L, :US, or :LS")
        return new{R,Vector{R},type}(dim, Vector{R}(undef, packedsize(dim)))
    end
end

const PackedMatrixUpper{R,V} = Union{PackedMatrix{R,V,:U},PackedMatrix{R,V,:US}}
const PackedMatrixLower{R,V} = Union{PackedMatrix{R,V,:L},PackedMatrix{R,V,:LS}}
const PackedMatrixUnscaled{R,V} = Union{PackedMatrix{R,V,:U},PackedMatrix{R,V,:L}}
const PackedMatrixScaled{R,V} = Union{PackedMatrix{R,V,:US},PackedMatrix{R,V,:LS}}

@inline packedsize(dim) = dim * (dim +1) ÷ 2
@inline rowcol_to_vec(P::PackedMatrixUpper, row, col) =
    (@boundscheck(1 ≤ row ≤ col || throw(BoundsError(P, (row, col)))); return col * (col -1) ÷ 2 + row)
@inline rowcol_to_vec(P::PackedMatrixLower, row, col) =
    (@boundscheck(1 ≤ col ≤ row || throw(BoundsError(P, (row, col)))); return (2P.dim - col) * (col -1) ÷ 2 + row)

packed_format(::PackedMatrix{R,V,Fmt}) where {R,V,Fmt} = Fmt
packed_isupper(::PackedMatrixUpper) = true
packed_isupper(::PackedMatrixLower) = false
packed_isupper(s::Symbol) = s == :U || s == :US
packed_islower(::PackedMatrixUpper) = false
packed_islower(::PackedMatrixLower) = true
packed_islower(s::Symbol) = s == :L || s == :LS
packed_isscaled(::PackedMatrixUnscaled) = false
packed_isscaled(::PackedMatrixScaled) = true
packed_isscaled(s::Symbol) = s == :US || s == :LS
packed_ulchar(x) = packed_isupper(x) ? 'U' : 'L'

struct PackedDiagonalIterator{Fmt}
    dim::Int
    k::UInt
    PackedDiagonalIterator(P::PackedMatrixUpper, k) = new{:U}(P.dim, abs(k))
    PackedDiagonalIterator(P::PackedMatrixLower, k) = new{:L}(P.dim, abs(k))
end
function Base.iterate(iter::PackedDiagonalIterator{:U})
    if iter.k ≥ iter.dim
        return nothing
    else
        j = iter.k * (iter.k +1) ÷ 2 +1
        return j, (j, iter.k +2)
    end
end
function Base.iterate(iter::PackedDiagonalIterator{:L})
    if iter.k ≥ iter.dim
        return nothing
    else
        j = iter.k +1
        return j, (j, iter.dim)
    end
end
function Base.iterate(iter::PackedDiagonalIterator{:U}, state)
    j, δ = state
    j += δ
    if δ > iter.dim
        return nothing
    else
        return j, (j, δ +1)
    end
end
function Base.iterate(::PackedDiagonalIterator{:L}, state)
    j, δ = state
    if isone(δ)
        return nothing
    else
        j += δ
        return j, (j, δ -1)
    end
end
Base.IteratorSize(::Type{PackedDiagonalIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{PackedDiagonalIterator}) = Base.HasEltype()
Base.eltype(::PackedDiagonalIterator) = Int
Base.length(iter::PackedDiagonalIterator) = iter.dim - iter.k
LinearAlgebra.diagind(A::PackedMatrix, k::Integer=0) = collect(PackedDiagonalIterator(A, k))
function LinearAlgebra.diag(A::PackedMatrix{R}, k::Integer=0) where {R}
    iter = PackedDiagonalIterator(A, k)
    diagonal = Vector{R}(undef, length(iter))
    for (i, idx) in enumerate(iter)
        @inbounds diagonal[i] = A[idx]
    end
    return diagonal
end
function LinearAlgebra.tr(A::PackedMatrix{R}) where {R}
    trace = zero(R)
    for idx in PackedDiagonalIterator(A, 0)
        @inbounds trace += A[idx]
    end
    return trace
end

function rmul_diags!(M::PackedMatrix{R}, factor::R) where {R}
    data = M.data
    for i in PackedDiagonalIterator(M, 0)
        @inbounds data[i] = data[i] * factor
    end
    return M
end
function rmul_offdiags!(M::PackedMatrix{R}, factor::R) where {R}
    diags = PackedDiagonalIterator(M, 0)
    data = M.data
    for (d₁, d₂) in zip(diags, Iterators.drop(diags, 1))
        @inbounds rmul!(@view(data[d₁+1:d₂-1]), factor)
    end
    return M
end
function lmul_diags!(M::PackedMatrix{R}, factor::R) where {R}
    data = M.data
    for i in PackedDiagonalIterator(M, 0)
        @inbounds data[i] = factor * data[i]
    end
    return M
end
function lmul_offdiags!(M::PackedMatrix{R}, factor::R) where {R}
    diags = PackedDiagonalIterator(M, 0)
    data = M.data
    for (d₁, d₂) in zip(diags, Iterators.drop(diags, 1))
        @inbounds lmul!(@view(data[d₁+1:d₂-1]), factor)
    end
    return M
end

packed_scale!(P::PackedMatrixScaled) = P
packed_scale!(P::PackedMatrixUnscaled{R}) where {R} =
    rmul_offdiags!(PackedMatrix(P.dim, P.data, packed_isupper(P) ? :US : :LS), sqrt(R(2)))
packed_unscale!(P::PackedMatrixScaled{R}) where {R} =
    rmul_offdiags!(PackedMatrix(P.dim, P.data, packed_isupper(P) ? :U : :L), sqrt(inv(R(2))))
packed_unscale!(P::PackedMatrixUnscaled) = P

function PackedMatrix(P::PackedMatrix, format::Symbol=packed_format(P))
    packed_isupper(P) == packed_isupper(format) || error("Changing the storage direction is not supported at the moment")
    if packed_isscaled(format)
        return packed_scale!(copy(P))
    else
        return packed_unscale!(copy(P))
    end
end

Base.size(P::PackedMatrix) = (P.dim, P.dim)
LinearAlgebra.checksquare(P::PackedMatrix) = P.dim
Base.eltype(::PackedMatrix{R}) where {R} = R
Base.@propagate_inbounds Base.getindex(P::PackedMatrix, idx) = P.data[idx]
Base.@propagate_inbounds Base.getindex(P::PackedMatrix{R,V,:U}, row, col) where {R,V} =
    P.data[@inbounds rowcol_to_vec(P, min(row, col), max(row, col))]
Base.@propagate_inbounds function Base.getindex(P::PackedMatrix{R,V,:US}, row, col) where {R,V}
    val = P.data[@inbounds rowcol_to_vec(P, min(row, col), max(row, col))]
    return row == col ? val : sqrt(inv(R(2))) * val
end
Base.@propagate_inbounds Base.getindex(P::PackedMatrix{R,V,:L}, row, col) where {R,V} =
    P.data[@inbounds rowcol_to_vec(P, max(row, col), min(row, col))]
Base.@propagate_inbounds function Base.getindex(P::PackedMatrix{R,V,:LS}, row, col) where {R,V}
    val = P.data[@inbounds rowcol_to_vec(P, max(row, col), min(row, col))]
    return row == col ? val : sqrt(inv(R(2))) * val
end
Base.@propagate_inbounds Base.setindex!(P::PackedMatrix, X, idx::Union{Integer,LinearIndices}) = P.data[idx] = X
Base.@propagate_inbounds Base.setindex!(P::PackedMatrix{R,V,:U}, X, row, col) where {R,V} =
    P.data[@inbounds rowcol_to_vec(P, min(row, col), max(row, col))] = X
Base.@propagate_inbounds function Base.setindex!(P::PackedMatrix{R,V,:US}, X, row, col) where {R,V}
    if row != col
        X *= sqrt(R(2))
    end
    P.data[@inbounds rowcol_to_vec(P, min(row, col), max(row, col))] = X
end
Base.@propagate_inbounds Base.setindex!(P::PackedMatrix{R,V,:L}, X, row, col) where {R,V} =
    P.data[@inbounds rowcol_to_vec(P, max(row, col), min(row, col))] = X
Base.@propagate_inbounds function Base.setindex!(P::PackedMatrix{R,V,:LS}, X, row, col) where {R,V}
    if row != col
        X *= sqrt(R(2))
    end
    P.data[@inbounds rowcol_to_vec(P, max(row, col), min(row, col))] = X
end
Base.IndexStyle(::PackedMatrix) = IndexLinear()
Base.IteratorSize(::Type{<:PackedMatrix}) = Base.HasLength()
Base.iterate(P::PackedMatrix, args...) = iterate(P.data, args...)
Base.length(P::PackedMatrix) = length(P.data)
Base.collect(P::PackedMatrix) = collect(P.data)
Base.fill!(P::PackedMatrix{R}, x::R) where {R} = fill!(P.data, x)
Base.copy(P::PackedMatrix) = PackedMatrix(P.dim, copy(P.data), packed_format(P))
for cp in (:copy!, :copyto!)
    for Fmt in (:(:U), :(:L), :(:US), :(:LS)) # we need to be so specific, as there would be ambiguities with the copyto!s below
        @eval begin
            Base.@propagate_inbounds function Base.$cp(dst::PackedMatrix{R1,V1,$Fmt}, src::PackedMatrix{R2,V2,$Fmt}) where {R1,V1,R2,V2}
                $cp(dst.data, src.data)
                return dst
            end
        end
    end
    @eval begin
        Base.@propagate_inbounds function Base.$cp(dst::PackedMatrix{R}, src::AbstractVector{R}) where {R}
            $cp(dst.data, src)
            return dst
        end
        Base.@propagate_inbounds Base.$cp(dst::AbstractVector{R}, src::PackedMatrix{R}) where {R} = $cp(dst, src.data)
    end
end
function Base.copyto!(dst::PackedMatrix{R,V,:U}, src::AbstractMatrix{R}) where {R,V}
    j = 1
    for (i, col) in enumerate(eachcol(src))
        @inbounds copyto!(dst.data, j, col, 1, i)
        j += i
    end
    return dst
end
function Base.copyto!(dst::PackedMatrix{R,V,:US}, src::AbstractMatrix{R}) where {R,V}
    j = 1
    @inbounds for (i, col) in enumerate(eachcol(src))
        @views dst.data[j:j+i-2] .= sqrt(R(2)) .* col[1:i-1]
        dst.data[j+i-1] = col[i]
        j += i
    end
    return dst
end
function Base.copyto!(dst::PackedMatrix{R,V,:L}, src::AbstractMatrix{R}) where {R,V}
    j = 1
    l = src.dim
    for (i, col) in enumerate(eachcol(src))
        @inbounds copyto!(dst.data, j, col, i, l)
        j += l
        l -= 1
    end
end
function Base.copyto!(dst::PackedMatrix{R,V,:LS}, src::AbstractMatrix{R}) where {R,V}
    j = 1
    l = src.dim
    @inbounds for (i, col) in enumerate(eachcol(src))
        dst.data[j] = col[i]
        @views dst.data[j+1:j+l-1] .= sqrt(R(2)) .* col[i+1:i+l-1]
        j += l
        l -= 1
    end
end
function Base.copy!(dst::PackedMatrix{R}, src::AbstractMatrix{R}) where {R}
    @boundscheck checkbounds(src, axes(dst)...)
    return copyto!(dst, src)
end
function Base.similar(P::PackedMatrix{R}, ::Type{T}=eltype(P), dims::NTuple{2,Int}=size(P); format::Symbol=packed_format(P)) where {R,T}
    ==(dims...) || error("Packed matrices must be square")
    dim = first(dims)
    return PackedMatrix(dim, similar(P.data, T, packedsize(dim)), format)
end
LinearAlgebra.vec(P::PackedMatrix) = P.data
Base.convert(T::Type{<:Ptr}, P::PackedMatrix) = convert(T, P.data)
Base.unsafe_convert(T::Type{Ptr{R}}, P::PackedMatrix{R}) where {R} = Base.unsafe_convert(T, P.data)
Base.reshape(P::PackedMatrix, ::Val{1}) = P.data # controversial? But it allows taking views with linear indices appropriately.

# Broadcasting
# idea: allow the PackedMatrix to work either as a vector or as a matrix, depending on the broadcasting context.
# However, this doesn't seem to work in all cases, so instead of providing a partially buggy implementation, we default to
# always treating the PackedMatrix as equivalent to its vector in broadcasting.
# We still perform a check so that we don't combine incompatible formats.
Base.broadcastable(P::PackedMatrix) = PackedMatrixBroadcasting{typeof(P)}(P)
struct PackedMatrixBroadcasting{PM<:PackedMatrix} <: AbstractVector{eltype(PM)}
    data::PM
end
Base.size(PB::PackedMatrixBroadcasting) = size(PB.data.data)
Base.axes(PB::PackedMatrixBroadcasting) = axes(PB.data.data)
Base.getindex(PB::PackedMatrixBroadcasting, ind::Int) = PB.data.data[ind]
Base.setindex!(PB::PackedMatrixBroadcasting, val, ind::Int) = PB.data.data[ind] = val
#=@inline Broadcast.combine_axes(A::PackedMatrixBroadcasting, B) = try
    return Broadcast.broadcast_shape(axes(A), axes(B))
catch
    return Broadcast.broadcast_shape(axes(A.data), axes(B))
end
@inline Broadcast.combine_axes(A, B::PackedMatrixBroadcasting) = try
    return Broadcast.broadcast_shape(axes(A), axes(B))
catch
    return Broadcast.broadcast_shape(axes(A), axes(B.data))
end
@inline Broadcast.combine_axes(A::PackedMatrixBroadcasting, B::PackedMatrixBroadcasting) = Broadcast.broadcast_shape(axes(A), axes(B))=#

struct PackedMatrixStyle{Fmt} <: Broadcast.AbstractArrayStyle{1} end
PackedMatrixStyle{Fmt}(::Val{1}) where {Fmt} = PackedMatrixStyle{Fmt}()
PackedMatrixStyle{Fmt}(::Val{2}) where {Fmt} = error("Broadcasting a PackedMatrix will only work on the vectorized data")
#PackedMatrixStyle{Fmt}(::Val{2}) where {Fmt} = PackedMatrixGenericStyle{Fmt}()
Base.BroadcastStyle(::Type{<:Union{<:PackedMatrixBroadcasting{<:PackedMatrix{R,V,Fmt} where {R,V}},<:PackedMatrix{R,V,Fmt} where {R,V}}}) where {Fmt} =
    PackedMatrixStyle{Fmt}()
Base.similar(bc::Broadcast.Broadcasted{P}, ::Type{T}) where {T,Fmt,P<:PackedMatrixStyle{Fmt}} =
    PackedMatrix{T}(undef, (isqrt(1 + 8length(bc)) -1) ÷ 2, Fmt)
    #=similar(find_pm(bc).data, T)
find_pm(bc::Base.Broadcast.Broadcasted) = find_pm(bc.args)
find_pm(args::Tuple) = find_pm(find_pm(args[1]), Base.tail(args))
find_pm(E::Broadcast.Extruded{<:PackedMatrixBroadcasting}) = E.x
find_pm(x) = x
find_pm(::Tuple{}) = nothing
find_pm(P::PackedMatrixBroadcasting, ::Any) = P
find_pm(::Any, rest) = find_pm(rest)=#
@inline function Base.copyto!(dest::PackedMatrix{R,V,Fmt} where {R,V}, bc::Broadcast.Broadcasted{PackedMatrixStyle{Fmt}}) where {Fmt}
    axes(dest.data) == axes(bc) || Broadcast.throwdm(axes(dest.data), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = Broadcast.preprocess(dest, bc)
    # Performance may vary depending on whether `@inbounds` is placed outside the
    # for loop or not. (cf. https://github.com/JuliaLang/julia/issues/38086)
    @inbounds @simd for I in eachindex(bc′)
        dest[I] = bc′[I]
    end
    return dest
end
Base.Broadcast.materialize!(s::PackedMatrixStyle{Fmt}, dest::PackedMatrix{R,V,Fmt} where {R,V}, bc::Broadcast.Broadcasted{PackedMatrixStyle{Fmt}}) where {Fmt} =
    (Base.Broadcast.materialize!(s, dest.data, bc); return dest)
Base.BroadcastStyle(::PackedMatrixStyle{Fmt}, ::PackedMatrixStyle{Fmt}) where {Fmt} = PackedMatrixStyle{Fmt}()
Base.BroadcastStyle(::PackedMatrixStyle, ::PackedMatrixStyle) =
    error("Packed matrices with different formats cannot be combined")
#=struct PackedMatrixGenericStyle{Fmt} <: Broadcast.AbstractArrayStyle{2} end
PackedMatrixGenericStyle{Fmt}(::Val{1}) where {Fmt} = PackedMatrixStyle{Fmt}()
PackedMatrixGenericStyle{Fmt}(::Val{2}) where {Fmt} = PackedMatrixGenericStyle{Fmt}()
Base.BroadcastStyle(::PackedMatrixStyle{Fmt}, ::Broadcast.DefaultArrayStyle{2}) where {Fmt} = PackedMatrixGenericStyle{Fmt}()
Base.BroadcastStyle(::PackedMatrixGenericStyle{Fmt}, ::Broadcast.DefaultArrayStyle{1}) where {Fmt} = PackedMatrixStyle{Fmt}()
Base.similar(bc::Broadcast.Broadcasted{<:PackedMatrixGenericStyle}, ::Type{T}) where {T} = similar(Array{T}, axes(bc))
@inline function Base.copyto!(dest::PackedMatrix{R,V,Fmt} where {R,V}, bc::Broadcast.Broadcasted{PackedMatrixGenericStyle{Fmt}}) where {Fmt}
    axes(dest) == axes(bc) || Broadcast.throwdm(axes(dest), axes(bc))
    if bc.f === identity && bc.args isa Tuple{AbstractArray}
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = Broadcast.preprocess(dest, bc)
    @inbounds @simd for I in eachindex(IndexCartesian(), bc′) # linear indices will be out of bounds!
        dest[I] = bc′[I]
    end
    return dest
end=#

chkpacked(n::Integer, AP::AbstractVector) = 2length(AP) == n * (n +1) || error("Packed storage length does not match dimension")
for (spev, spevd, spevx, pptrf, spmv, spr, tpttr, trttp, gemmt, elty) in
    ((:dspev_,:dspevd_,:dspevx_,:dpptrf_,:dspmv_,:dspr_,:dtpttr_,:dtrttp_,:dgemmt_,:Float64),
     (:sspev_,:sspevd_,:sspevx_,:spptrf_,:sspmv_,:sspr_,:stpttr_,:strttp_,:sgemmt_,:Float32))
     @eval begin
        #       SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, UPLO
        #       INTEGER            INFO, LDZ, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
        function spev!(jobz::AbstractChar, uplo::AbstractChar, n::Integer, AP::AbstractVector{$elty})
            BLAS.chkuplo(uplo)
            LinearAlgebra.chkstride1(AP)
            chkpacked(n, AP)
            W     = similar(AP, $elty, n)
            work  = Vector{$elty}(undef, 3n)
            info  = Ref{BLAS.BlasInt}()
            ldz   = n
            if jobz == 'N'
                Z = similar(AP, $elty, ldz, 0)
            else
                Z = similar(AP, $elty, ldz, n)
            end
            ccall((BLAS.@blasfunc($spev), BLAS.libblastrampoline), Cvoid,
                  (Ref{UInt8},        # JOBZ
                   Ref{UInt8},        # UPLO
                   Ref{BLAS.BlasInt}, # N
                   Ptr{$elty},        # AP
                   Ptr{$elty},        # W
                   Ptr{$elty},        # Z
                   Ref{BLAS.BlasInt}, # LDZ
                   Ptr{$elty},        # WORK
                   Ptr{BLAS.BlasInt}, # INFO
                   Clong,             # length(JOBZ)
                   Clong),            # length(UPLO)
                jobz, uplo, n, AP, W, Z, ldz, work, info, 1, 1)
            LAPACK.chklapackerror(info[])
            jobz == 'V' ? (W, Z) : W
        end

        #       SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
        #      $                   ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, RANGE, UPLO
        #       INTEGER            IL, INFO, IU, LDZ, M, N
        #       DOUBLE PRECISION   ABSTOL, VL, VU
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IFAIL( * ), IWORK( * )
        #       DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
        function spevx!(jobz::AbstractChar, range::AbstractChar, uplo::AbstractChar, n::Integer, AP::AbstractVector{$elty},
                        vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)
            BLAS.chkuplo(uplo)
            LinearAlgebra.chkstride1(AP)
            chkpacked(n, AP)
            if range == 'I' && !(1 <= il <= iu <= n)
                throw(ArgumentError("illegal choice of eigenvalue indices (il = $il, iu = $iu), which must be between 1 and n = $n"))
            end
            if range == 'V' && vl >= vu
                throw(ArgumentError("lower boundary, $vl, must be less than upper boundary, $vu"))
            end
            m = Ref{BLAS.BlasInt}()
            W = similar(AP, $elty, n)
            ldz = n
            if jobz == 'N'
                Z = similar(AP, $elty, ldz, 0)
            elseif jobz == 'V'
                Z = similar(AP, $elty, ldz, n)
            end
            work   = Vector{$elty}(undef, 8n)
            iwork  = Vector{BLAS.BlasInt}(undef, 5n)
            ifail  = Vector{BLAS.BlasInt}(undef, n)
            info   = Ref{BLAS.BlasInt}()
            ccall((BLAS.@blasfunc($spevx), BLAS.libblastrampoline), Cvoid,
                  (Ref{UInt8},        # JOBZ
                   Ref{UInt8},        # RANGE
                   Ref{UInt8},        # UPLO
                   Ref{BLAS.BlasInt}, # N
                   Ptr{$elty},        # AP
                   Ref{$elty},        # VL
                   Ref{$elty},        # VU
                   Ref{BLAS.BlasInt}, # IL
                   Ref{BLAS.BlasInt}, # IU
                   Ref{$elty},        # ABSTOL
                   Ptr{BLAS.BlasInt}, # M
                   Ptr{$elty},        # W
                   Ptr{$elty},        # Z
                   Ref{BLAS.BlasInt}, # LDZ
                   Ptr{$elty},        # WORK
                   Ptr{BLAS.BlasInt}, # IWORK
                   Ptr{BLAS.BlasInt}, # IFAIL
                   Ptr{BLAS.BlasInt}, # INFO
                   Clong,             # length(JOBZ)
                   Clong,             # length(RANGE)
                   Clong),            # length(UPLO)
                jobz, range, uplo, n, AP, vl, vu, il, iu, abstol, m, W, Z, max(1, ldz), work, iwork, ifail, info, 1, 1, 1)
            LAPACK.chklapackerror(info[])
            W[1:m[]], Z[:,1:(jobz == 'V' ? m[] : 0)]
        end
        spevx!(jobz::AbstractChar, n::Integer, A::AbstractVector{$elty}) =
            spevx!(jobz, 'A', 'U', n, A, 0.0, 0.0, 0, 0, -1.0)

        #           SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK,
        #      $                   LWORK, IWORK, LIWORK, INFO )
        # *     .. Scalar Arguments ..
        # CHARACTER          JOBZ, UPLO
        # INTEGER            INFO, LDZ, LIWORK, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        # INTEGER            IWORK( * )
        # DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
        function spevd!(jobz::AbstractChar, uplo::AbstractChar, n::Integer, AP::AbstractVector{$elty})
            BLAS.chkuplo(uplo)
            LinearAlgebra.chkstride1(AP)
            chkpacked(n, AP)
            W      = similar(AP, $elty, n)
            ldz    = n
            if jobz == 'N'
                Z = similar(AP, $elty, ldz, 0)
            else
                Z = similar(AP, $elty, ldz, n)
            end
            work   = Vector{$elty}(undef, 1)
            lwork  = BLAS.BlasInt(-1)
            iwork  = Vector{LinearAlgebra.BlasInt}(undef, 1)
            liwork = BLAS.BlasInt(-1)
            info   = Ref{BLAS.BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((BLAS.@blasfunc($spevd), BLAS.libblastrampoline), Cvoid,
                      (Ref{UInt8},        # JOBZ
                       Ref{UInt8},        # UPLO
                       Ref{BLAS.BlasInt}, # N
                       Ptr{$elty},        # AP
                       Ptr{$elty},        # W
                       Ptr{$elty},        # Z
                       Ref{BLAS.BlasInt}, # LDZ
                       Ptr{$elty},        # WORK
                       Ref{BLAS.BlasInt}, # LWORK
                       Ptr{BLAS.BlasInt}, # IWORK
                       Ref{BLAS.BlasInt}, # LIWORK
                       Ptr{BLAS.BlasInt}, # INFO
                       Clong,             # length(JOBZ)
                       Clong),            # length(UPLO)
                      jobz, uplo, n, AP, W, Z, ldz, work, lwork, iwork, liwork, info, 1, 1)
                LAPACK.chklapackerror(info[])
                if i == 1
                    lwork = BLAS.BlasInt(work[1])
                    resize!(work, lwork)
                    liwork = BLAS.BlasInt(iwork[1])
                    resize!(iwork, liwork)
                end
            end
            jobz == 'V' ? (W, Z) : W
        end

        # SUBROUTINE DPPTRF( UPLO, N, AP, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   AP( * )
        function pptrf!(uplo::AbstractChar, n::Integer, AP::AbstractVector{$elty})
            LinearAlgebra.require_one_based_indexing(AP)
            LinearAlgebra.chkstride1(AP)
            BLAS.chkuplo(uplo)
            chkpacked(n, AP)
            info = Ref{BLAS.BlasInt}()
            ccall((BLAS.@blasfunc($pptrf), BLAS.libblastrampoline), Cvoid,
                  (Ref{UInt8},        # UPLO
                   Ref{BLAS.BlasInt}, # N
                   Ptr{$elty},        # AP
                   Ptr{BLAS.BlasInt}, # INFO
                   Clong),            # length(UPLO)
                  uplo, n, AP, info, 1)
            LAPACK.chkargsok(info[])
            #info[] > 0 means the leading minor of order info[] is not positive definite
            #ordinarily, throw Exception here, but return error code here
            #this simplifies isposdef! and factorize
            return AP, info[] # info stored in Cholesky
        end

        #           SUBROUTINE DSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
        # *     .. Scalar Arguments ..
        #         DOUBLE PRECISION ALPHA,BETA
        #         INTEGER INCX,INCY,N
        #         CHARACTER UPLO
        #   *     ..
        #   *     .. Array Arguments ..
        #         DOUBLE PRECISION AP(*),X(*),Y(*)
        function spmv!(uplo::AbstractChar, α::Union{$elty,Bool}, AP::AbstractVector{$elty}, x::AbstractVector{$elty},
                       β::Union{$elty,Bool}, y::AbstractVector{$elty})
            BLAS.chkuplo(uplo)
            LinearAlgebra.chkstride1(AP)
            LinearAlgebra.require_one_based_indexing(AP, x, y)
            n = length(x)
            length(y) == n || ArgumentError("Input and output dimensions must match")
            chkpacked(n, AP)
            px, stx = BLAS.vec_pointer_stride(x, ArgumentError("input vector with 0 stride is not allowed"))
            py, sty = BLAS.vec_pointer_stride(y, ArgumentError("dest vector with 0 stride is not allowed"))
            GC.@preserve x AP ccall((BLAS.@blasfunc($spmv), BLAS.libblastrampoline), Cvoid,
                (Ref{UInt8},        # UPLO
                 Ref{BLAS.BlasInt}, # N
                 Ref{$elty},        # ALPHA
                 Ptr{$elty},        # AP
                 Ptr{$elty},        # X
                 Ref{BLAS.BlasInt}, # INCX
                 Ref{$elty},        # BETA
                 Ptr{$elty},        # Y
                 Ref{BLAS.BlasInt}, # INCY
                 Clong),            # length(UPLO)
                uplo, n, α, AP, px, stx, β, py, sty, 1)
            return y
        end

        #           SUBROUTINE DSPR( UPLO, N, ALPHA, X, INCX, AP )
        # *     .. Scalar Arguments ..
        #       DOUBLE PRECISION ALPHA
        #       INTEGER INCX,N
        #       CHARACTER UPLO
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION AP(*),X(*)
        function spr!(uplo::AbstractChar, α::$elty, x::AbstractVector{$elty}, AP::AbstractVector{$elty})
            BLAS.chkuplo(uplo)
            LinearAlgebra.chkstride1(AP)
            LinearAlgebra.require_one_based_indexing(AP, x)
            n = length(x)
            chkpacked(n, AP)
            px, stx = BLAS.vec_pointer_stride(x, ArgumentError("input vector with 0 stride is not allowed"))
            GC.@preserve x ccall((BLAS.@blasfunc($spr), BLAS.libblastrampoline), Cvoid,
                (Ref{UInt8},        # UPLO
                 Ref{BLAS.BlasInt}, # N
                 Ref{$elty},        # ALPHA
                 Ptr{$elty},        # X
                 Ref{BLAS.BlasInt}, # INCX
                 Ptr{$elty},        # AP
                 Clong),            # length(UPLO)
                uplo, n, α, px, stx, AP, 1)
            return AP
        end

        #           SUBROUTINE DTPTTR( UPLO, N, AP, A, LDA, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, N, LDA
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), AP( * )
        function tpttr!(uplo::AbstractChar, AP::AbstractVector{$elty}, A::StridedArray{$elty})
            BLAS.chkuplo(uplo)
            LinearAlgebra.chkstride1(AP)
            LinearAlgebra.require_one_based_indexing(AP, A)
            n = LinearAlgebra.checksquare(A)
            chkpacked(n, AP)
            lda = Base.stride(A, 2)
            info = Ref{BLAS.BlasInt}()
            ccall((BLAS.@blasfunc($tpttr), BLAS.libblastrampoline), Cvoid,
                (Ref{UInt8},        # UPLO
                 Ref{BLAS.BlasInt}, # N
                 Ptr{$elty},        # AP
                 Ptr{$elty},        # A
                 Ref{BLAS.BlasInt}, # LDA
                 Ptr{BLAS.BlasInt}, # INFO
                 Clong),            # length(UPLO)
                uplo, n, AP, A, lda, info, 1)
            LAPACK.chklapackerror(info[])
            return A
        end

        #           SUBROUTINE DTRTTP( UPLO, N, A, LDA, AP, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, N, LDA
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), AP( * )
        function trttp!(uplo::AbstractChar, A::StridedArray{$elty}, AP::AbstractVector{$elty})
            BLAS.chkuplo(uplo)
            LinearAlgebra.chkstride1(AP)
            LinearAlgebra.require_one_based_indexing(AP, A)
            n = LinearAlgebra.checksquare(A)
            chkpacked(n, AP)
            lda = Base.stride(A, 2)
            info = Ref{BLAS.BlasInt}()
            ccall((BLAS.@blasfunc($trttp), BLAS.libblastrampoline), Cvoid,
                (Ref{UInt8},        # UPLO
                 Ref{BLAS.BlasInt}, # N
                 Ptr{$elty},        # A
                 Ref{BLAS.BlasInt}, # LDA
                 Ptr{$elty},        # AP
                 Ptr{BLAS.BlasInt}, # INFO
                 Clong),            # length(UPLO)
                uplo, n, A, lda, AP, info, 1)
            LAPACK.chklapackerror(info[])
            return AP
        end

        # Not part of the reference BLAS (yet), but present in many implementations, also Julia 1.9
        # Requires at least OpenBlas 0.3.22, but Julia currently maximally ships with 0.3.17. But with this construction, we
        # make the function available anyway - just load MKL before this package. When the function is not found, the more
        # inefficient gemm! is called.
        #           SUBROUTINE DGEMMT( UPLO, TRANSA, TRANSB, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
        #       .. Scalar Arguments ..
        #       DOUBLE PRECISION ALPHA,BETA
        #       INTEGER K,LDA,LDB,LDC,N
        #       CHARACTER TRANSA,TRANSB, UPLO
        #       ..
        #       .. Array Arguments ..
        #       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
        if BLAS.lbt_get_forward(BLAS.@blasfunc($gemmt),
            BLAS.USE_BLAS64 ? BLAS.LBT_INTERFACE_ILP64 : BLAS.LBT_INTERFACE_LP64) ∈ (BLAS.lbt_get_default_func(), Ptr{Cvoid}(-1))
            function gemmt!(uplo::AbstractChar, transA::AbstractChar, transB::AbstractChar, alpha::Union{($elty), Bool},
                A::AbstractVecOrMat{$elty}, B::AbstractVecOrMat{$elty}, beta::Union{($elty), Bool},
                C::AbstractVecOrMat{$elty})
                return BLAS.gemm!(transA, transB, alpha, A, B, beta, C)
            end
        else
            function gemmt!(uplo::AbstractChar, transA::AbstractChar, transB::AbstractChar, alpha::Union{($elty), Bool},
                            A::AbstractVecOrMat{$elty}, B::AbstractVecOrMat{$elty}, beta::Union{($elty), Bool},
                            C::AbstractVecOrMat{$elty})
                LinearAlgebra.require_one_based_indexing(A, B, C)
                BLAS.chkuplo(uplo)
                m = size(A, transA == 'N' ? 1 : 2)
                ka = size(A, transA == 'N' ? 2 : 1)
                kb = size(B, transB == 'N' ? 1 : 2)
                n = size(B, transB == 'N' ? 2 : 1)
                if ka != kb || m != LinearAlgebra.checksquare(C) || n != m
                    throw(DimensionMismatch(lazy"A has size ($m,$ka), B has size ($kb,$n), C has size $(size(C))"))
                end
                LinearAlgebra.chkstride1(A)
                LinearAlgebra.chkstride1(B)
                LinearAlgebra.chkstride1(C)
                ccall((BLAS.@blasfunc($gemmt), BLAS.libblastrampoline), Cvoid,
                    (Ref{UInt8},        # UPLO
                    Ref{UInt8},        # TRANSA
                    Ref{UInt8},        # TRANSB
                    Ref{BLAS.BlasInt}, # N
                    Ref{BLAS.BlasInt}, # K
                    Ref{$elty},        # ALPHA
                    Ptr{$elty},        # A
                    Ref{BLAS.BlasInt}, # LDA
                    Ptr{$elty},         # B
                    Ref{BLAS.BlasInt}, # LDB
                    Ref{$elty},        # BETA
                    Ptr{$elty},        # C
                    Ref{BLAS.BlasInt}, # LDC
                    Clong,             # length(UPLO)
                    Clong,             # length(TRANSA)
                    Clong,             # length(TRANSB)
                    ),
                    uplo, transA, transB, n, ka, alpha, A, max(1, Base.stride(A, 2)), B, max(1, Base.stride(B,2)), beta, C,
                    max(1, Base.stride(C, 2)), 1, 1, 1)
                C
            end
        end
    end
end

function tpttr!(uplo::AbstractChar, AP::AbstractVector{T}, A::AbstractMatrix{T}) where {T}
    BLAS.chkuplo(uplo)
    n = LinearAlgebra.checksquare(A)
    chkpacked(n, AP)
    if uplo == 'L' || uplo == 'l'
        k = 0
        @inbounds for j in 1:n, i in j:n
            k += 1
            A[i, j] = AP[k]
        end
    else
        k = 0
        @inbounds for j in 1:n, i in 1:j
            k += 1
            A[i, j] = AP[k]
        end
    end
    return A
end

function trttp!(uplo::AbstractChar, A::AbstractMatrix{T}, AP::AbstractVector{T}) where {T}
    BLAS.chkuplo(uplo)
    n = LinearAlgebra.checksquare(A)
    chkpacked(n, AP)
    if uplo == 'L' || uplo == 'l'
        k = 0
        @inbounds for j in 1:n, i in j:n
            k += 1
            AP[k] = A[i, j]
        end
    else
        k = 0
        @inbounds for j in 1:n, i in 1:j
            k += 1
            AP[k] = A[i, j]
        end
    end
    return AP
end

LinearAlgebra.axpy!(a::Number, x::PackedMatrix, y::AbstractVector) = LinearAlgebra.axpy!(a, x.data, y)
LinearAlgebra.axpy!(a::Number, x::AbstractVector, y::PackedMatrix) = LinearAlgebra.axpy!(a, x, y.data)
LinearAlgebra.axpy!(a::Number, x::PackedMatrix, y::PackedMatrix) = LinearAlgebra.axpy!(a, x.data, y.data)

# These are for multiplications where the matrix is directly interpreted as a vector.
LinearAlgebra.mul!(C::PackedMatrix, A::AbstractMatrix, B::AbstractVector, α::Number, β::Number) = mul!(C.data, A, B, α, β)
LinearAlgebra.mul!(C::PackedMatrix, A::AbstractMatrix, B::PackedMatrix, α::Number, β::Number) = mul!(C.data, A, B.data, α, β)
LinearAlgebra.mul!(C::AbstractVector, A::AbstractMatrix, B::PackedMatrix, α::Number, β::Number) = mul!(C, A, B.data, α, β)
# But we need to be very careful: if A above was a PackedMatrix, they all transform to Vector, PackedMatrix, Vector, which will
# (potentially) call the mul! below that performs symmetric multiplication - which is not what should be done.
LinearAlgebra.mul!(::Union{<:AbstractVector,PackedMatrix}, ::PackedMatrix, ::Union{<:AbstractVector,<:PackedMatrix}, ::Number,
    ::Number) = throw(DimensionMismatch())
# And this is the symmetric multiplication. BLAS only offers it for the packed format. So either we write a helper that
# - allocates a copy and unscales it (bad - mul! shouldn't allocate) or
# - unscales the matrix, multiplies and rescales (bad - floating point not necessarily reversible; not threadsafe)
# or we leave this unsupported.
LinearAlgebra.mul!(C::AbstractVector, A::PackedMatrixUnscaled, B::AbstractVector, α::Number, β::Number) =
    spmv!(packed_ulchar(A), α, A.data, B, β, C)
# Also a wrapper for rank-one updates. Here, we are allowed to mutate AP, so we un-scale it. But beware, the type may change!
function spr!(α, x::AbstractVector, AP::PackedMatrix)
    APu = packed_unscale!(AP)
    spr!(packed_ulchar(APu), α, x, APu.data)
    return APu
end
function Base.Matrix{R}(A::PackedMatrixUnscaled{R}) where {R}
    result = Matrix{R}(undef, A.dim, A.dim)
    tpttr!(packed_ulchar(A), A.data, result)
    return Symmetric(result, packed_isupper(A) ? :U : :L)
end
function Base.Matrix{R}(A::PackedMatrix{R,V,:US}) where {R,V}
    result = Matrix{R}(undef, A.dim, A.dim)
    tpttr!('U', A.data, result)
    for j in 2:A.dim
        @inbounds rmul!(@view(result[1:j-1, j]), sqrt(inv(R(2))))
    end
    return Symmetric(result, :U)
end
function Base.Matrix{R}(A::PackedMatrix{R,V,:LS}) where {R,V}
    result = Matrix{R}(undef, A.dim, A.dim)
    tpttr!('L', A.data, result)
    for j in 1:A.dim-1
        @inbounds rmul!(@view(result[j+1:end, j]), sqrt(inv(R(2))))
    end
    return Symmetric(result, :L)
end
function PackedMatrix(A::Symmetric{R,<:AbstractMatrix{R}}) where {R}
    result = PackedMatrix{R}(undef, size(A, 1), A.uplo == 'U' ? :U : :L)
    trttp!(A.uplo, parent(A), result.data)
    return result
end

function LinearAlgebra.dot(A::PackedMatrix{R1,V,Fmt} where {V}, B::PackedMatrix{R2,V,Fmt} where {V}) where {R1,R2,Fmt}
    A.dim == B.dim || error("Matrices must have same dimensions")
    if packed_isscaled(Fmt)
        return dot(A.data, B.data)
    else
        result = zero(promote_type(R1, R2))
        diags = PackedDiagonalIterator(A, 0)
        cur_diag = iterate(diags)
        i = 1
        @inbounds while !isnothing(cur_diag)
            @simd for j in i:cur_diag[1]-1
                result += 2conj(A[j]) * B[j]
            end
            result += conj(A[cur_diag[1]]) * B[cur_diag[1]]
            i = cur_diag[1] +1
            cur_diag = iterate(diags, cur_diag[2])
        end
        return result
    end
end
function LinearAlgebra.dot(A::PackedMatrix{R1,<:SparseVector,Fmt}, B::PackedMatrix{R2,V,Fmt} where {V}) where {R1,R2,Fmt}
    A.dim == B.dim || error("Matrices must have same dimensions")
    if packed_isscaled(Fmt)
        return dot(A.data, B.data)
    else
        result = zero(promote_type(R1, R2))
        nzs = rowvals(A.data)
        vs = nonzeros(A.data)
        diags = PackedDiagonalIterator(A, 0)
        cur_diag = iterate(diags)
        isnothing(cur_diag) && return result
        @inbounds for i in 1:length(nzs) # for @simd, cannot iterate over (j, v)
            j = nzs[i]
            v = vs[i]
            while cur_diag[1] < j
                cur_diag = iterate(diags, cur_diag[2])
                if isnothing(cur_diag)
                    cur_diag = (typemax(Int), cur_diag[2])
                end
            end
            if cur_diag[1] == j
                result += conj(v) * B[j]
            else
                result += 2conj(v) * B[j]
            end
        end
        return result
    end
end
LinearAlgebra.dot(A::PackedMatrix{R}, B::PackedMatrix{R,<:SparseVector}) where {R} = conj(dot(B, A))
function LinearAlgebra.dot(A::PackedMatrix{R1,<:SparseVector,Fmt}, B::PackedMatrix{R2,<:SparseVector,Fmt}) where {R1,R2,Fmt}
    A.dim == B.dim || error("Matrices must have same dimensions")
    if packed_isscaled(Fmt)
        return dot(A.data, B.data)
    else
        result = zero(promote_type(R1, R2))
        nzA = rowvals(A.data)
        nzB = rowvals(B.data)
        vA = nonzeros(A.data)
        vB = nonzeros(B.data)
        iAmax = length(nzA)
        iBmax = length(nzB)
        iA = 1
        iB = 1
        diags = PackedDiagonalIterator(A, 0)
        cur_diag = iterate(diags)
        isnothing(cur_diag) && return result
        @inbounds while iA ≤ iAmax && iB ≤ iBmax
            pA = nzA[iA]
            pB = nzB[iB]
            while pA > pB
                iB += 1
                iB > iBmax && return result
                pB = nzB[iB]
            end
            while pB > pA
                iA += 1
                iA > iAmax && return result
                pA = nzA[iA]
            end
            @assert pA == pB
            while cur_diag[1] < pA
                cur_diag = iterate(diags, cur_diag[2])
                if isnothing(cur_diag)
                    cur_diag = (typemax(Int), cur_diag[2])
                end
            end
            if cur_diag[1] == pA
                result += conj(vA[iA]) * vB[iB]
            else
                result += 2conj(vA[iA]) * vB[iB]
            end
            iA += 1
            iB += 1
        end
        return result
    end
end

function normapply(f, pm::PackedMatrix{R}, init::T=zero(R)) where {R,T}
    result::T = init
    diags = PackedDiagonalIterator(pm, 0)
    cur_diag = iterate(diags)
    i = 1
    @inbounds while !isnothing(cur_diag)
        @simd for j in i:cur_diag[1]-1
            if packed_isscaled(pm)
                result = f(result, pm.data[j] * sqrt(inv(R(2))), false)
            else
                result = f(result, pm.data[j], false)
            end
        end
        result = f(result, pm.data[cur_diag[1]], true)
        i = cur_diag[1] +1
        cur_diag = iterate(diags, cur_diag[2])
    end
    return result
end
function normapply(f, pm::PackedMatrix{R,<:SparseVector}) where {R}
    nzs = rowvals(pm.data)
    vs = nonzeros(pm.data)
    diags = PackedDiagonalIterator(pm, 0)
    cur_diag = iterate(diags)
    isnothing(cur_diag) && return zero(R)
    result = zero(R)
    @inbounds for i in 1:length(nzs)
        j = nzs[i]
        while cur_diag[1] < j
            cur_diag = iterate(diags, cur_diag[2])
            if isnothing(cur_diag)
                cur_diag = (typemax(Int), cur_diag[2])
            end
        end
        if cur_diag[1] == j
            result = f(result, vs[i], true)
        elseif packed_isscaled(pm)
            result = f(result, vs[i] * sqrt(inv(R(2))), false)
        else
            result = f(result, vs[i], false)
        end
    end
    return result
end
LinearAlgebra.norm2(pm::PackedMatrixUnscaled) = sqrt(normapply((Σ, x, diag) -> Σ + (diag ? abs(x)^2 : 2abs(x)^2), pm))
LinearAlgebra.norm2(pm::PackedMatrixScaled) = LinearAlgebra.norm2(pm.data)
LinearAlgebra.norm2(pm::PackedMatrixScaled{R,<:SparseVector}) where {R} = LinearAlgebra.norm2(nonzeros(pm.data))
LinearAlgebra.norm1(pm::PackedMatrix) = normapply((Σ, x, diag) -> Σ + (diag ? abs(x) : 2abs(x)), pm)
LinearAlgebra.normInf(pm::PackedMatrixUnscaled) = LinearAlgebra.normInf(pm.data)
LinearAlgebra.normInf(pm::PackedMatrixScaled) = normapply((m, x, diag) -> max(m, abs(x)), pm)
Base._simple_count(f, pm::PackedMatrix, init::T) where {T} =
    normapply((Σ, x, diag) -> f(x) ? (diag ? Σ + one(T) : Σ + one(T) + one(T)) : Σ, pm, init)
LinearAlgebra.normMinusInf(pm::PackedMatrixUnscaled) = LinearAlgebra.normMinusInf(pm.data)
LinearAlgebra.normMinusInf(pm::PackedMatrixScaled{R}) where {R} = normapply((m, x, diag) -> min(m, abs(x)), pm, R(Inf))
LinearAlgebra.normp(pm::PackedMatrix, p) = normapply((Σ, x, diag) -> Σ + (diag ? abs(x)^p : 2abs(x)^p), pm)^(1/p)

LinearAlgebra.eigen!(pm::PackedMatrixUnscaled{R}) where {R<:Real} =
    Eigen(spevd!('V', packed_ulchar(pm), pm.dim, pm.data)...)
function LinearAlgebra.eigen!(pm::PackedMatrixScaled{R}) where {R<:Real}
    fac = sqrt(R(2))
    eval, evec = spevd!('V', packed_ulchar(pm), pm.dim, rmul_diags!(pm, fac).data)
    return Eigen(rmul!(eval, inv(fac)), evec)
end
LinearAlgebra.eigvals(pm::PackedMatrixUnscaled{R}) where {R<:Real} =
    spevd!('N', packed_ulchar(pm), pm.dim, copy(pm.data))
function LinearAlgebra.eigvals(pm::PackedMatrixScaled{R}) where {R<:Real}
    fac = sqrt(R(2))
    return rmul!(spevd!('N', packed_ulchar(pm), pm.dim, rmul_diags!(copy(pm), fac).data), inv(fac))
end
LinearAlgebra.eigen!(pm::PackedMatrixUnscaled{R}, vl::R, vu::R) where {R<:Real} =
    Eigen(spevx!('V', 'V', packed_ulchar(pm), pm.dim, pm.data, vl, vu, 0, 0, -one(R))...)
function LinearAlgebra.eigen!(pm::PackedMatrixScaled{R}, vl::R, vu::R) where {R<:Real}
    fac = sqrt(R(2))
    eval, evec = spevx!('V', 'V', packed_ulchar(pm), pm.dim, rmul_diags!(pm, fac).data, vl * fac, vu * fac, 0, 0, -one(R))
    return Eigen(rmul!(eval, inv(fac)), evec)
end
LinearAlgebra.eigen!(pm::PackedMatrixUnscaled{R}, range::UnitRange) where {R<:Real} =
    Eigen(spevx!('V', 'I', packed_ulchar(pm), pm.dim, pm.data, zero(R), zero(R), range.start, range.stop, -one(R))...)
function LinearAlgebra.eigen!(pm::PackedMatrixScaled{R}, range::UnitRange) where {R<:Real}
    fac = sqrt(R(2))
    eval, evec = spevx!('V', 'I', packed_ulchar(pm), pm.dim, rmul_diags!(pm, fac).data, zero(R), zero(R), range.start,
        range.stop, -one(R))
    return Eigen(rmul!(eval, inv(fac)), evec)
end
LinearAlgebra.eigvals!(pm::PackedMatrixUnscaled{R}, vl::R, vu::R) where {R<:Real} =
    spevx!('N', 'V', packed_ulchar(pm), pm.dim, pm.data, vl, vu, 0, 0, -one(R))[1]
function LinearAlgebra.eigvals!(pm::PackedMatrixScaled{R}, vl::R, vu::R) where {R<:Real}
    fac = sqrt(R(2))
    return rmul!(spevx!('N', 'V', packed_ulchar(pm), pm.dim, rmul_diags!(pm, fac).data, vl * fac, vu * fac, 0, 0, -one(R))[1],
        inv(fac))
end
LinearAlgebra.eigvals!(pm::PackedMatrixUnscaled{R}, range::UnitRange) where {R<:Real} =
    spevx!('N', 'I', packed_ulchar(pm), pm.dim, pm.data, zero(R), zero(R), range.start, range.stop, -one(R))[1]
function LinearAlgebra.eigvals!(pm::PackedMatrixScaled{R}, range::UnitRange) where {R<:Real}
    fac = sqrt(R(2))
    return rmul!(spevx!('N', 'I', packed_ulchar(pm), pm.dim, rmul_diags!(pm, fac).data, zero(R), zero(R), range.start,
        range.stop, -one(R))[1], inv(fac))
end
eigmin!(pm::PackedMatrix) = eigvals!(pm, 1:1)[1]
eigmax!(pm::PackedMatrix) = eigvals!(pm, pm.dim:pm.dim)[1]
LinearAlgebra.eigmin(pm::PackedMatrix) = eigmin!(copy(pm))
LinearAlgebra.eigmax(pm::PackedMatrix) = eigmax!(copy(pm))
function LinearAlgebra.cholesky!(pm::PackedMatrix{R}, ::NoPivot = NoPivot(); shift::R = zero(R), check::Bool = true) where {R<:Real}
    if !iszero(shift)
        for i in PackedDiagonalIterator(pm, 0)
            @inbounds pm[i] += shift
        end
    end
    C, info = pptrf!(packed_ulchar(pm), pm.dim, packed_unscale!(pm).data)
    check && LinearAlgebra.checkpositivedefinite(info)
    return Cholesky(PackedMatrix(pm.dim, C, packed_isupper(pm) ? :U : :L), packed_ulchar(pm), info)
end
LinearAlgebra.isposdef(pm::PackedMatrix{R}, tol::R=zero(R)) where {R<:Real} =
    isposdef(cholesky!(copy(pm), shift=tol, check=false))

end
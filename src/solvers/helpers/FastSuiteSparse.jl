# This file contains a couple of enhancements for the Julia interface to SuiteSparse by using functions that allow for passing
# preallocated data.
import SuiteSparse: CHOLMOD, SPQR

struct EfficientCholmod{T,F<:Factorization{T}} <: Factorization{T}
    cholmod::F
    Y::Ref{Ptr{CHOLMOD.cholmod_dense}}
    E::Ref{Ptr{CHOLMOD.cholmod_dense}}

    function EfficientCholmod(f::CHOLMOD.Factor)
        result = new{eltype(f),typeof(f)}(f, Ref{Ptr{CHOLMOD.cholmod_dense}}(C_NULL), Ref{Ptr{CHOLMOD.cholmod_dense}}(C_NULL))
        finalizer(f) do _
            let result = result
                common = CHOLMOD.getcommon()
                CHOLMOD.cholmod_free_dense(result.Y, common)
                CHOLMOD.cholmod_free_dense(result.E, common)
            end
        end
        return result
    end
end

function cholmod_dense_wrap(A::StridedVecOrMat{T}) where {T}
    out = CHOLMOD.cholmod_dense()
    out.nrow = size(A, 1)
    out.ncol = size(A, 2)
    out.d = Base.stride(A, 2)
    out.nzmax = out.d * out.ncol
    out.x = pointer(A)
    out.z = C_NULL
    out.xtype = CHOLMOD.xtyp(T)
    out.dtype = CHOLMOD.getcommon()[].dtype
    return out
end

function LinearAlgebra.ldiv!(out::StridedVecOrMat{T}, lhs::EfficientCholmod{T,<:Factorization{T}},
    rhs::StridedVecOrMat{T}) where {T<:CHOLMOD.VTypes}
    F = lhs.cholmod
    if size(F,1) != size(rhs,1)
        throw(DimensionMismatch("LHS and RHS should have the same number of rows. " *
            "LHS has $(size(F,1)) rows, but RHS has $(size(rhs,1)) rows."))
    end
    if size(F,2) != size(out,1) # F is square, but anyway...
        throw(DimensionMismatch("Output must have the same number of rows as LHS has columns. " *
            "LHS has $(size(F,2)) columns, but output has $(size(out,1)) rows."))
    end
    if size(out,2) != size(rhs,2)
        throw(DimensionMismatch("Output and RHS should have the same number of column." *
            "Output has $(size(out,2)) columns, but RHS has $(size(rhs,2)) columns."))
    end
    if !CHOLMOD.issuccess(F)
        s = unsafe_load(pointer(F))
        if s.is_ll == 1
            throw(LinearAlgebra.PosDefException(s.minor))
        else
            throw(LinearAlgebra.ZeroPivotException(s.minor))
        end
    end
    # Casting via CHOLMOD.Dense{T}(rhs) will actually copy around all the data.
    B = cholmod_dense_wrap(rhs)
    X = cholmod_dense_wrap(out)
    pX = Ref(X)
    GC.@preserve pX CHOLMOD.cholmod_l_solve2(
        CHOLMOD.CHOLMOD_A, # system to solve
        F,                 # factorization to use
        Ref(B), C_NULL,    # right-hand-side, dense and sparse
        Ref(Base.unsafe_convert(Ptr{CHOLMOD.cholmod_dense}, pX)),
        C_NULL,            # solution, dense and sparse
        lhs.Y, lhs.E,      # workspace
        CHOLMOD.getcommon()
    ) == 0 && error("Error in cholmod")
    return out
end

LinearAlgebra.ldiv!(lhs::EfficientCholmod{T,<:Factorization{T}}, rhs::StridedVecOrMat{T}) where {T<:CHOLMOD.VTypes} =
    ldiv!(rhs, lhs, rhs) # this appears to work, although it is not documented

# for unknown reasons, there is no zero-allocation sparse-QR ldiv!, although it would be most easy to implement. We first
# implement the interface to the vector permutation function LAPMR in BLAS.
for (lapmr, elty) in ((:dlapmr_,:Float64), (:slapmr_,:Float32))
    @eval begin
       #       SUBROUTINE DLAPMR( FORWRD, M, N, X, LDX, K )
       # *     .. Scalar Arguments ..
       #       LOGICAL            FORWRD
       #       INTEGER            LDX, M, N
       # *     ..
       # *     .. Array Arguments ..
       #       INTEGER            K( * )
       #       DOUBLE PRECISION   X( LDX, * )
       # *     ..
       function lapmr!(forward::Bool, X::StridedMatrix{$elty}, K::AbstractVector{Int})
            LinearAlgebra.chkstride1(K)
            if length(K) != size(X, 1)
                throw(ArgumentError("permutation must have $(size(X, 1)) elements, found $(length(K))"))
            end
            ccall((BLAS.@blasfunc($lapmr), BLAS.libblastrampoline), Cvoid,
                 (Ref{BLAS.BlasInt}, # FORWRD
                  Ref{BLAS.BlasInt}, # M
                  Ref{BLAS.BlasInt}, # N
                  Ptr{$elty},        # X
                  Ref{BLAS.BlasInt}, # LDX
                  Ptr{BLAS.BlasInt}),# K
               forward, size(X, 1), size(X, 2), X, Base.stride(X, 2), K)
            return X
       end
   end
end

function LinearAlgebra.ldiv!(X::StridedVecOrMat{T}, F::SPQR.QRSparse, B::StridedVecOrMat{T}) where {T}
    if size(F, 1) != size(B, 1)
        throw(DimensionMismatch("size(F) = $(size(F)) but size(B) = $(size(B))"))
    end
    if size(X, 1) != max(size(F, 2), size(B, 1))
        throw(DimensionMismatch("$(size(X, 1)) output rows but need $(max(size(F, 2), size(B, 1)))"))
    end
    if size(X, 2) != size(B, 2)
        throw(DimensionMismatch("size(X) = $(size(X)) but size(B) = $(size(B))"))
    end

    # The rank of F equal might be reduced
    rnk = rank(F)

    # Fill will zeros. These will eventually become the zeros in the basic solution
    # fill!(X, 0)
    # Apply left permutation to the solution and store in X
    for j in 1:size(B, 2)
        for i in 1:length(F.rpivinv)
            @inbounds X[F.rpivinv[i], j] = B[i, j]
        end
    end

    # Make a view into X corresponding to the size of B
    X0 = view(X, 1:size(B, 1), :)

    # Apply Q' to B
    LinearAlgebra.lmul!(adjoint(F.Q), X0)

    # Zero out to get basic solution
    X[rnk + 1:end, :] .= 0

    # Solve R*X = B
    LinearAlgebra.ldiv!(UpperTriangular(F.R[Base.OneTo(rnk), Base.OneTo(rnk)]),
                        view(X0, Base.OneTo(rnk), :))

    # Apply right permutation and extract solution from X
    # NB: cpiv == [] if SPQR was called with ORDERING_FIXED
    Xout = @view(X[1:size(F,2), :])
    if !isempty(F.cpiv)
        if X isa AbstractMatrix
            lapmr!(false, Xout, F.cpiv)
        else
            invpermute!(Xout, F.cpiv)
        end
    end
    return Xout
end

function LinearAlgebra.ldiv!(F::SPQR.QRSparse, B::StridedVecOrMat{T}) where {T}
    if size(F, 1) != size(B, 1)
        throw(DimensionMismatch("size(F) = $(size(F)) but size(B) = $(size(B))"))
    end
    if size(B, 1) < size(F, 2)
        throw(DimensionMismatch("$(size(B, 1)) output rows but need $(size(F, 2))"))
    end

    # The rank of F equal might be reduced
    rnk = rank(F)

    # Fill will zeros. These will eventually become the zeros in the basic solution
    # fill!(X, 0)
    # Apply left permutation to the solution and store in X
    if !isempty(F.rpivinv)
        if B isa AbstractMatrix
            lapmr!(false, B, F.rpivinv)
        else
            invpermute!(B, F.rpivinv)
        end
    end

    # Apply Q' to B
    LinearAlgebra.lmul!(adjoint(F.Q), B)

    # Zero out to get basic solution
    B[rnk + 1:end, :] .= 0

    # Solve R*X = B
    LinearAlgebra.ldiv!(UpperTriangular(F.R[Base.OneTo(rnk), Base.OneTo(rnk)]),
                        view(B, Base.OneTo(rnk), :))

    # Apply right permutation and extract solution from X
    # NB: cpiv == [] if SPQR was called with ORDERING_FIXED
    if !isempty(F.cpiv)
        if B isa AbstractMatrix
            lapmr!(false, B, F.cpiv)
        else
            invpermute!(B, F.cpiv)
        end
    end
    return B
end
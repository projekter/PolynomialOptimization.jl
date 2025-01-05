# Here, we provide a function to delete columns from a matrix.
# Julia doesn't make this easy, as a Matrix is managed by the C core. Therefore, we hack ourselves access into the internals
# (which still doesn't work well unless the Matrix was constructed by unsafe_wrap of a malloc'ed memory).
# But starting from version 1.11, this changed and all vector-like stuff is much more Julia-managed. We can therefore make
# things much simpler and neater.
if VERSION < v"1.11"
    # mutable struct jl_array_t
    #     data::Ptr{Cvoid}
    #     length::UInt
    #     flags::UInt16 # how:2, ndims:9, pooled:1, ptrarray:1, hasptr:1, isshared:1, isaligned:1
    #     elsize::UInt16
    #     offset::UInt32
    #     nrows::UInt
    #     maxsize::UInt
    # end
    const jl_array_flags_isaligned = 0b1_0_0_0_0_000000000_00
    const jl_array_flags_isshared  = 0b0_1_0_0_0_000000000_00
    const jl_array_flags_hasptr    = 0b0_0_1_0_0_000000000_00
    const jl_array_flags_ptrarray  = 0b0_0_0_1_0_000000000_00
    const jl_array_flags_pooled    = 0b0_0_0_0_1_000000000_00
    const jl_array_flags_ndims     = 0b0_0_0_0_0_111111111_00
    const jl_array_flags_how       = 0b0_0_0_0_0_000000000_11

    function matrix_delete_end!(a::Matrix, ncols::Integer)
        iszero(ncols) && return a
        s = size(a)
        1 ≤ ncols ≤ s[2] || throw(BoundsError(a, (1, ncols)))
        ap = pointer_from_objref(a)
        flags = unsafe_load(Ptr{UInt16}(ap) + sizeof(Ptr) + sizeof(UInt))
        if flags & jl_array_flags_how != 2
            return a[:, 1:s[2]-ncols]
            # This is bad, but there's not really a better way. If we were to simulate _array_del_end, this would just change
            # the internally stored size and not free any memory. This is due to the fact that 1d arrays have the maxsize field
            # that holds the actual size, which is used for overallocation. However, for 2d, maxsize doesn't exist and instead
            # must hold the number of rows. There is no overallocation. So we don't have any other choice but resizing the
            # memory, which could be done by calling jl_array_shrink, only that it is not available. Instead, it is wrapped
            # into jl_array_sizehint, which does a check whether we actually decrease by a noteworthy amount, else it just
            # returns. Hence, it might not do the job. But even if it did the job, it would just allocate a new buffer and copy
            # the data over - so we can just do it in a canonical way.
        else
            # Here, things are different. malloc-allocated pointers are deferred to realloc in jl_array_shrink. We still cannot
            # call this function, but we can reproduce its behavior.
            l = length(a)
            newl = (s[2] - ncols) * s[1]
            elsz = sizeof(eltype(a))
            # jl_array_isbitsunion: !flags.ptrarray && jl_is_uniontype(eltype)
            isbitsunion = iszero(flags & jl_array_flags_ptrarray) && Base.isbitsunion(eltype(a))
            if isbitsunion
                newbytes = newl * (elsz +1)
                oldnbytes = l * (elsz +1)
                typetagdata = Libc.malloc(newl)
                ccall(:memcpy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
                    typetagdata,
                    ccall(:jl_array_typetagdata, Ptr{Cchar}, (Any,), a),
                    newl
                )
            else
                newbytes = newl * elsz
                oldnbytes = l * elsz
                if elsz == 1
                    newbytes += 1
                    oldnbytes += 1
                end
            end
            data_p = Ptr{Ptr{Cvoid}}(ap)
            length_p = Ptr{UInt}(ap) + sizeof(Ptr)
            offset_p = Ptr{UInt32}(length_p) + sizeof(UInt) + 2sizeof(UInt16)
            ncols_p = Ptr{UInt}(offset_p) + sizeof(UInt32) + sizeof(UInt)

            oldoffsnb = unsafe_load(offset_p) * elsz
            originalptr = unsafe_load(data_p) - oldoffsnb
            # change data
            unsafe_store!(
                data_p,
                ccall(:jl_gc_managed_realloc, Ptr{Cvoid}, (Ptr{Cvoid}, Csize_t, Csize_t, Cint, Any),
                    originalptr, newbytes, oldnbytes, flags & jl_array_flags_isaligned, a) + oldoffsnb
            )
            # change length
            unsafe_store!(length_p, newl)
            # change number of columns
            unsafe_store!(ncols_p, s[2] - ncols)
            if isbitsunion
                newtypetagdata = ccall(:jl_array_typetagdata, Ptr{Cchar}, (Any,), a)
                ccall(:memcpy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), newtypetagdata, typetagdata, newl)
                Libc.free(typetagdata)
            end
            return a
        end
    end

    function resizable_array(T::Type, dims...)
        a = unsafe_wrap(Array, Ptr{T}(Libc.malloc(*(dims...) * sizeof(T))), dims, own=true)
        # by default, a is shared already, which is clearly not correct here
        ap = pointer_from_objref(a)
        flagsptr = Ptr{UInt16}(ap) + sizeof(Ptr) + sizeof(UInt)
        unsafe_store!(flagsptr, unsafe_load(flagsptr) & ~jl_array_flags_isshared)
        return a
    end
    resizable_copy(a::Array{T}) where {T} = @inbounds copyto!(resizable_array(T, size(a)...), a)
else
    # TODO: check
    function matrix_delete_end!(a::Matrix, ncols::Integer)
        iszero(ncols) && return a
        s = size(a)
        1 ≤ ncols ≤ s[2] || throw(BoundsError(a, (1, ncols)))
        l = length(a)
        newl = (s[2] - ncols) * s[1]
        ref = a.ref
        for i in newl+1:l
            @inbounds Base._unsetindex!(GenericMemoryRef(ref, i))
        end
        setfield!(a, :size, (s[1], (s[2] - ncols)))
        # We also want to free the memory, so we copy a behavior similar to sizehint. However, shrinking involves copying, so
        # only do it if it's worth doing.
        if l - newl > l ÷ 8
            mem = ref.mem
            newmem = Base.array_new_memory(mem, newl)
            newref = GenericMemoryRef(newmem)
            unsafe_copyto!(newref, ref, newl)
            setfield(a, :ref, newref)
        end
        return a
    end

    resizable_array(T::Type, dims...) = Array{T}(undef, dims)
    resizable_copy(a::Array) = copy(a)
end

"""
    keepcol!(A::Union{<:Matrix,<:SparseMatrixCSC}, m::AbstractVector{Bool})

Only keep the columns in the matrix `A` that are `true` in the vector `m` (which must have the same length as `size(A, 2)`).
The output is logically the same as `A[:, m]`; however, note that the data underlying `A` is mutated: this function does not
create a copy. The function returns the new matrix; the old one should no longer be used (it might become invalid, as the
sparse matrix struct is immutable).
"""
function keepcol! end

function keepcol!(A::SparseMatrixCSC, m::AbstractVector{Bool})
    d = size(A, 2)
    length(m) == d || throw(BoundsError(A, (1, m)))
    colΔ = 0
    colptr, rowval, nzval = SparseArrays.getcolptr(A), rowvals(A), nonzeros(A)
    i = 1
    last = 1
    from = 1
    @inbounds while from ≤ d && !m[from]
        from += 1
    end
    to = from
    @inbounds while to ≤ d
        if m[to]
            to += 1
        else
            len = to - from
            datalen = colptr[to] - colptr[from]
            colΔ += colptr[from] - colptr[last]
            if !iszero(colΔ)
                copyto!(rowval, colptr[from] - colΔ, rowval, colptr[from], datalen)
                copyto!(nzval, colptr[from] - colΔ, nzval, colptr[from], datalen)
            end
            isone(from) || (colptr[i:i+len-1] .= colptr[from:to-1] .- colΔ)
            i += len
            last = to
            from = to +1
            while from ≤ d && !m[from]
                from += 1
            end
            to = from
        end
    end
    @inbounds if from ≤ d
        len = to - from
        datalen = colptr[to] - colptr[from]
        colΔ += colptr[from] - colptr[last]
        if !iszero(colΔ)
            copyto!(rowval, colptr[from] - colΔ, rowval, colptr[from], datalen)
            copyto!(nzval, colptr[from] - colΔ, nzval, colptr[from], datalen)
        end
        isone(from) || (colptr[i:i+len-1] .= colptr[from:to-1] .- colΔ)
        i += len
    end
    colptr[i] = colptr[end] - colΔ
    # this just removes the entries from the array
    resize!
    Base._deleteend!(rowval, colΔ)
    Base._deleteend!(nzval, colΔ)
    Base._deleteend!(colptr, d +1 - i)
    # but here, we actually make it smaller (if this makes sense)
    sizehint!(rowval, length(rowval))
    sizehint!(nzval, length(nzval))
    sizehint!(colptr, length(colptr))
    return SparseMatrixCSC{eltype(nzval),eltype(colptr)}(size(A, 1), length(colptr) -1, colptr, rowval, nzval)
end

function keepcol!(A::Matrix, m::AbstractVector{Bool})
    inc, d = size(A)
    length(m) == d || throw(BoundsError(A, m))
    i = 1
    from = 1
    @inbounds while from ≤ d && !m[from]
        from += 1
    end
    to = from
    @inbounds while to ≤ d
        if m[to]
            to += 1
        else
            len = to - from # [from, to)
            isone(from) || copyto!(A, (i -1) * inc +1, A, (from -1) * inc +1, len * inc)
            i += len
            from = to +1
            while from ≤ d && !m[from]
                from += 1
            end
            to = from
        end
    end
    if from ≤ d
        len = to - from # [from, to), as now to = d +1
        @inbounds isone(from) || copyto!(A, (i -1) * inc +1, A, (from -1) * inc +1, len * inc)
        i += len
    end
    return matrix_delete_end!(A, d - i +1)
end
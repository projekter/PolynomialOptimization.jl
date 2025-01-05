export Indvals, IndvalsIterator, PSDMatrixCartesian, PSDIndextypeMatrixCartesian, PSDIndextypeVector, psd_indextype

"""
    Indvals{T,V}

Supertype for an iterable that returns a `Tuple{T,V}` on iteration, where the first is a variable/constraint index and the
second its coefficient in the constraint matrix. The parameters of the type correspond to those in [`AbstractSolver`](@ref).
The properties `indices` and `values` can be accessed and will give `AbstractVector`s of the appropriate type. Note that the
fields should only be used if an iterative approach is not feasible, as they might be constructed on-demand (this will only
happen for the first two indvals in the standard quadratic cone, all other elements can be accessed with zero cost).
"""
struct Indvals{T,V,Z}
    z::Z

    function Indvals(indices::AbstractVector{T}, values::AbstractVector{V}) where {T,V<:Real}
        @assert(length(indices) == length(values))
        z = zip(indices, values)
        new{T,V,typeof(z)}(z)
    end
end

function Base.getproperty(iv::Indvals, f::Symbol)
    if f === :indices
        return getfield(iv, :z).is[1]
    elseif f === :values
        return getfield(iv, :z).is[2]
    else
        return getfield(iv, f)
    end
end

Base.IteratorSize(::Type{<:Indvals}) = Base.HasLength()
Base.IteratorEltype(::Type{<:Indvals}) = Base.HasEltype()
Base.eltype(::Type{Indvals{T,V}}) where {T,V} = Tuple{T,V}
Base.length(iv::Indvals) = length(iv.z)
Base.iterate(iv::Indvals, args...) = iterate(iv.z, args...)

"""
    IndvalsIterator{T,V}

An iterable that returns consecutive elements in a vectorized PSD cone.
This type stores a vector of indices and values together with information about the length of the individual subsequences.
Iterating through it will give [`Indvals`](@ref) that contain views into the indices and the values.
The vector of indices is available via `SparseArrays.rowvals`, the vector of values via `SparseArrays.nonzeros`, and the
lengths of the subsequences via `Base.index_lengths`.

See also [`PSDIndextypeVector`](@ref).
"""
struct IndvalsIterator{T,V,L,VT<:AbstractVector{T},VV<:AbstractVector{V}}
    indices::VT
    values::VV
    lens::L

    function IndvalsIterator(indices::AbstractVector{T}, values::AbstractVector{V}, len::I) where {T,V,I<:Integer}
        length(indices) == length(values) || throw(ArgumentError("Invalid IndvalsIterator construction"))
        new{T,V,I,typeof(indices),typeof(values)}(indices, values, len)
    end
    IndvalsIterator(::Unsafe, indices::AbstractVector{T}, values::AbstractVector{V}, lens::L) where {T,V,L<:Union{<:Integer,<:AbstractVector{<:Integer}}} =
        new{T,V,L,typeof(indices),typeof(values)}(indices, values, lens)
    function IndvalsIterator(indices::AbstractVector{T}, values::AbstractVector{V}, lens::L) where {T,V,L<:AbstractVector{<:Integer}}
        length(indices) == length(values) == sum(lens, init=0)|| error("Invalid IndvalsIterator construction")
        new{T,V,L,typeof(indices),typeof(values)}(indices, values, lens)
    end
end

Base.IteratorSize(::Type{<:IndvalsIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:IndvalsIterator}) = Base.HasEltype()
Base.eltype(::Type{IndvalsIterator{T,V,<:Any}}) where {T,V} =
    Tuple{SubArray{T,1,Vector{T},Tuple{UnitRange{Int}},true},SubArray{V,1,Vector{V},Tuple{UnitRange{Int}},true}}
Base.length(psdi::IndvalsIterator{<:Any,<:Any,<:Integer}) = length(psdi.indices) ÷ psdi.lens
Base.length(psdi::IndvalsIterator{<:Any,<:Any,<:AbstractVector{<:Integer}}) = length(psdi.lens)
@inline function Base.iterate(psdi::IndvalsIterator{<:Any,<:Any,<:Integer}, state=1)
    endpos = state + psdi.lens -1
    if endpos ≤ length(psdi.indices)
        @inbounds return Indvals(view(psdi.indices, state:endpos), view(psdi.values, state:endpos)), endpos +1
    else
        return nothing
    end
end
@inline function Base.iterate(psdi::IndvalsIterator{<:Any,<:Any,<:AbstractVector{<:Integer}}, state=(1, 1))
    state[2] ≤ length(psdi.lens) || return nothing
    startpos = state[1]
    endpos = @inbounds(startpos + psdi.lens[state[2]] -1)
    @inbounds return Indvals(view(psdi.indices, startpos:endpos), view(psdi.values, startpos:endpos)), (endpos +1, state[2] +1)
end
SparseArrays.rowvals(psdi::IndvalsIterator) = psdi.indices
SparseArrays.nonzeros(psdi::IndvalsIterator) = psdi.values
Base.index_lengths(psdi::IndvalsIterator{<:Any,<:Any,<:Integer}) = ConstantVector(psdi.lens, length(psdi))
Base.index_lengths(psdi::IndvalsIterator{<:Any,<:Any,<:AbstractVector{<:Integer}}) = psdi.lens

"""
    PSDMatrixCartesian

An iterable that returns matrix elements of a PSD cone.
Iterating through it will give `Pair`s with keys that contain the [`mindex`](@ref) of the monomial, and a 3-Tuple of
`AbstractVector`s as the values which contain the row and column indices together with their coefficients to describe where
the monomial appears. Note that the vectors will be reused in every iteration.

See also [`PSDIndextypeMatrixCartesian`](@ref).
"""
struct PSDMatrixCartesian{T,V,Tri,Offset}
    rowind::Vector{T}
    colind::Vector{T}
    nzval::Vector{V}
    trinumbers::Vector{T}
    symrows::FastVec{T}
    symcols::FastVec{T}

    function PSDMatrixCartesian{Offset}(dim::Int, tri::Symbol, rowind::Vector{T}, colind::Vector{T}, nzval::Vector{V}) where {T,V,Offset}
        tri === :L || tri === :U || error("Invalid triangle")
        length(rowind) == length(colind) == length(nzval) || error("Invalid PSDMatrixCartesian construction")
        sort_along!(colind, rowind, nzval)
        new{T,V,tri,
            if isbitstype(T) # be careful with large integers. But assume the offset will never exceed machine numbers
                T(Offset)
            elseif T <: Unsigned
                UInt(Offset)
            else
                Int(Offset)
            end
           }(rowind, colind, nzval, trisize.(zero(T):T(dim)-1), FastVec{T}(), FastVec{T}())
    end
end

Base.IteratorSize(::Type{<:PSDMatrixCartesian}) = Base.HasLength()
Base.IteratorEltype(::Type{<:PSDMatrixCartesian}) = Base.HasEltype()
Base.eltype(::Type{<:PSDMatrixCartesian{T,V}}) where {T,V} =
    Pair{T,Tuple{FastVec{T},FastVec{T},SubArray{V,1,Vector{V},Tuple{UnitRange{Int}},true}}}
function Base.length(psdi::PSDMatrixCartesian)
    isempty(psdi.colind) && return 0
    result = 1
    @inbounds for i in 2:length(psdi.colind)
        if psdi.colind[i] != psdi.colind[i-1]
            result += 1
        end
    end
    return result
end
@inline function Base.iterate(psdi::PSDMatrixCartesian{<:Any,<:Any,Tri,Offset}, (i, remaining)=(1, length(psdi.colind))) where {Tri,Offset}
    iszero(remaining) && return nothing
    rowind, colind, trinumbers, symrows, symcols = psdi.rowind, psdi.colind, psdi.trinumbers, psdi.symrows, psdi.symcols
    empty!(symrows)
    empty!(symcols)
    @inbounds curcol = colind[i]
    @inbounds while true
        idx = rowind[i]
        triindex = searchsortedlast(trinumbers, idx)
        # The input data is always in U format, apart from the fact that the nzval are already conjugated if Tri === :L
        push!(Tri === :L ? symrows : symcols, triindex - one(triindex) + Offset)
        push!(Tri === :L ? symcols : symrows, idx - trinumbers[triindex] + Offset)
        if isone(remaining) || colind[i+1] != curcol
            return Pair(curcol, (symrows, symcols, view(psdi.nzval, i-length(symrows)+1:i))), (i +1, remaining -1)
        end
        i += 1
        remaining -= 1
    end
end

"""
    PSDIndextypeMatrixCartesian(triangle, offset) <: PSDIndextype

The solver implements PSD matrix constraints by using a monolithic PSD matrix variable or an LMI-style representation.
Entries from the variable are obtained (or put into the LMI) by using a cartesian index of two integers of the type parameter
`T` of the [`AbstractSolver`](@ref). This index represents one triangle of the matrix (the lower if `triangle === :L`, the
upper if `triangle === :U`). The first entry has the index `(offset, offset)`, typically either `0` or `1`.

!!! info
    Note that while only one triangle is indexed, it is assumed that the solver will by default populate the other triangle in
    a completely symmetric way.
    This corresponds to the typical behavior of solvers that expose PSD variables and allow accessing elements in them via
    sparse symmetric matrices where only one triangle is given, but the other half is implicit.

See also [`PSDMatrixCartesian`](@ref).
"""
struct PSDIndextypeMatrixCartesian{Tri,Offset}
    function PSDIndextypeMatrixCartesian(triangle::Symbol, offset::Integer)
        triangle ∈ (:L, :U) || throw(MethodError(PSDIndextypeMatrixCartesian, (triangle, offset)))
        new{triangle,offset}()
    end
end

_get_offset(::PSDIndextypeMatrixCartesian{<:Any,Offset}) where {Offset} = Offset

"""
    PSDIndextypeVector(triangle) <: PSDIndextype

The solver implements PSD matrix constraints by demanding that the matrixization of a vector of decision variables be PSD.

If `triangle === :F`, the vector is formed by stacking all the columns of the matrix.

If `triangle === :L`, the columns of the lower triangle are assumed to be stacked _and scaled_, i.e., off-diagonal variables
that enter the cone are implicitly multiplied by ``1 / \\sqrt2`` in the matrix; so the coefficients will already be
premultiplied by ``\\sqrt2`` (for the [`add_constr_psd!`](@ref) case) or by ``1 / \\sqrt2`` (for the [`add_var_psd!`](@ref)
case).

If `triangle === :U`, the columns of the upper triangle are assumed to be stacked and scaled.

See also [`IndvalsIterator`](@ref).
"""
struct PSDIndextypeVector{Tri}
    function PSDIndextypeVector(triangle::Symbol)
        triangle ∈ (:L, :U, :F) || throw(MethodError(PSDIndextypeVector, (triangle,)))
        new{triangle}()
    end
end

"""
    PSDIndextype{Tri}

Union for all supported types in which a solver can represent a PSD matrix.

See also [`PSDIndextypeMatrixCartesian`](@ref), [`PSDIndextypeVector`](@ref).
"""
const PSDIndextype{Tri} = Union{<:PSDIndextypeMatrixCartesian{Tri},PSDIndextypeVector{Tri}}

"""
    psd_indextype(::AbstractSolver)

This function must indicate in which format the solver expects its data for PSD variables. The return type must be an instance
of a [`PSDIndextype`](@ref) subtype.

See also [`PSDIndextypeMatrixCartesian`](@ref), [`PSDIndextypeVector`](@ref).
"""
function psd_indextype end
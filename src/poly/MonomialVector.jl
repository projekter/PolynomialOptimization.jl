export SimpleMonomialVector, effective_nvariables, change_backend

struct SimpleMonomialVector{Nr,Nc,I<:Integer,E,T<:SimpleMonomial{Nr,Nc,I}} <: AbstractVector{T}
    data::E

    @inline function SimpleMonomialVector{Nr,Nc}(e::AbstractExponentsDegreeBounded{N,I}) where {Nr,Nc,N,I<:Integer}
        N == Nr + 2Nc || throw(MethodError(SimpleMonomialVector{Nr,Nc}, (e,)))
        index_counts(e, 0) # populate the cache
        new{Nr,Nc,I,typeof(e),SimpleMonomial{Nr,Nc,I,typeof(e)}}(e)
    end

    @inline function SimpleMonomialVector{Nr,Nc}(::Unsafe, e::AbstractExponents{N,I}, indices::AbstractVector{I}) where {Nr,Nc,N,I<:Integer}
        N == Nr + 2Nc || throw(MethodError(SimpleMonomialVector{Nr,Nc}, (unsafe, e, indices))) # compile-time check
        @inbounds degree_from_index(e, isempty(indices) ? one(I) : indices[end]) # populate the cache
        et = (e, indices)
        new{Nr,Nc,I,typeof(et),SimpleMonomial{Nr,Nc,I,typeof(e)}}(et)
    end
end

function SimpleMonomialVector{Nr,Nc}(e::AbstractExponents{<:Any,I}, indices::AbstractVector{I}) where {Nr,Nc,I<:Integer}
    _sortedallunique(sort!(indices)) || throw(ArgumentError("Monomial vector must not contain duplicates"))
    if !isempty(indices)
        @boundscheck @inbounds indices[begin] ≥ firstindex(e) || throw(BoundsError(e, indices[begin]))
        e isa AbstractExponentsDegreeBounded &&
            (@boundscheck @inbounds indices[end] ≤ lastindex(e) || throw(BoundsError(e, indices[end])))
    end
    return SimpleMonomialVector{Nr,Nc}(unsafe, e, indices)
end

SimpleMonomialVector{Nr,Nc}(::Unsafe, indices::AbstractVector{I}) where {Nr,Nc,I<:Integer} =
    SimpleMonomialVector{Nr,Nc}(unsafe, ExponentsAll{Nr+2Nc,I}(), indices)
SimpleMonomialVector{Nr,Nc}(indices::AbstractVector{I}) where {Nr,Nc,I<:Integer} =
    SimpleMonomialVector{Nr,Nc}(ExponentsAll{Nr+2Nc,I}(), indices)

"""
    SimpleMonomialVector{Nr,0[,I]}([e::AbstractExponents,]
        exponents_real::AbstractMatrix{<:Integer}, along...)
    SimpleMonomialVector{0,Nc[,I]}([e::AbstractExponents,]
        exponents_complex::AbstractMatrix{<:Integer},
        exponents_conj::AbstractMatrix{<:Integer}, along...)
    SimpleMonomialVector{Nr,Nc[,I]}([e::AbstractExponents,]
        exponents_real::AbstractMatrix{<:Integer},
        exponents_complex::AbstractMatrix{<:Integer},
        exponents_conj::AbstractMatrix{<:Integer}, along...)

Creates a monomial vector, where each column corresponds to one monomial and each row contains its exponents. The internal
representation will be made with respect to the exponent set `e`. If `e` is omitted, `ExponentsAll{Nr+2Nc,UInt}` is chosen by
default. Alternatively, all three methods may also be called with the index type `I` as a third type parameter, omitting `e`,
which then chooses `ExponentsAll{Nr+2Nc,I}` by default.
All matrices must have the same number of columns; complex and conjugate matrices must have the same number of rows.
The input will be sorted; if `along` are present, those vectors will be put in the same order as the inputs.
The input must not contain duplicates.
"""
function SimpleMonomialVector{Nr,0}(e::AbstractExponents{Nr,I}, exponents_real::AbstractMatrix{<:Integer}, along...) where {Nr,I<:Integer}
    size(exponents_real, 1) == Nr || throw(ArgumentError("Requested $Nr real variables, but got $(size(exponents_real, 1))"))
    exps = Vector{I}(undef, size(exponents_real, 2))
    for (i, col) in zip(eachindex(exps), eachcol(exponents_real))
        @inbounds exps[i] = exponents_to_index(e, col)
    end
    sort_along!(exps, along...)
    _sortedallunique(exps) || throw(ArgumentError("Monomial vector must not contain duplicates"))
    if !isempty(exps)
        @boundscheck @inbounds exps[begin] ≥ firstindex(e) || throw(BoundsError(e, exps[begin]))
        e isa AbstractExponentsDegreeBounded &&
            (@boundscheck @inbounds exps[end] ≤ lastindex(e) || throw(BoundsError(e, exps[end])))
    end
    return SimpleMonomialVector{Nr,0}(unsafe, e, exps)
end

function SimpleMonomialVector{0,Nc}(e::AbstractExponents{N,I}, exponents_complex::AbstractMatrix{<:Integer},
    exponents_conj::AbstractMatrix{<:Integer}, along...) where {Nc,N,I<:Integer}
    N == 2Nc || throw(MethodError(SimpleMonomialVector{0,Nc}, (e, exponents_complex, exponents_conj, along...)))
    size(exponents_complex, 1) == size(exponents_conj, 1) ||
        throw(ArgumentError("Complex and conjugate exponents lengths are different"))
    size(exponents_complex, 1) == Nc ||
        throw(ArgumentError("Requested $Nc complex variables, but got $(size(exponents_complex, 1))"))
    size(exponents_complex, 2) == size(exponents_conj, 2) ||
        throw(ArgumentError("Number of monomials is different"))
    exps = Vector{I}(undef, size(exponents_complex, 2))
    for (i, col) in zip(eachindex(exps), zip(eachcol(exponents_complex), eachcol(exponents_conj)))
        @inbounds exps[i] = exponents_to_index(e, (x[i] for i in 1:Nc for x in col))
    end
    sort_along!(exps, along...)
    _sortedallunique(exps) || throw(ArgumentError("Monomial vector must not contain duplicates"))
    if !isempty(exps)
        @boundscheck @inbounds exps[begin] ≥ firstindex(e) || throw(BoundsError(e, exps[begin]))
        e isa AbstractExponentsDegreeBounded &&
            (@boundscheck @inbounds exps[end] ≤ lastindex(e) || throw(BoundsError(e, exps[end])))
    end
    return SimpleMonomialVector{0,Nc}(unsafe, e, exps)
end

function SimpleMonomialVector{Nr,Nc}(e::AbstractExponents{N,I}, exponents_real::AbstractMatrix{<:Integer},
    exponents_complex::AbstractMatrix{<:Integer}, exponents_conj::AbstractMatrix{<:Integer}, along...) where {Nr,Nc,N,I<:Integer}
    N == Nr + 2Nc || throw(MethodError(SimpleMonomialVector{Nr,Nc}, (e, exponents_real, exponents_complex, exponents_conj,
                                                                     along...)))
    size(exponents_real, 1) == Nr || throw(ArgumentError("Requested $Nr real variables, but got $(size(exponents_real, 1))"))
    size(exponents_complex, 1) == size(exponents_conj, 1) ||
        throw(ArgumentError("Complex and conjugate exponents lengths are different"))
    size(exponents_complex, 1) == Nc ||
        throw(ArgumentError("Requested $Nc complex variables, but got $(size(exponents_complex, 1))"))
    size(exponents_real, 2) == size(exponents_complex, 2) == size(exponents_conj, 2) ||
        throw(ArgumentError("Number of monomials is different"))
    exps = Vector{I}(undef, size(exponents_real, 2))
    for (i, colreal, colcomplex, colconj) in zip(eachindex(exps), eachcol(exponents_real), eachcol(exponents_complex),
                                                    eachcol(exponents_conj))
        @inbounds exps[i] = exponents_to_index(e, OrderedExponents(colreal, colcomplex, colconj))
    end
    sort_along!(exps, along...)
    _sortedallunique(exps) || throw(ArgumentError("Monomial vector must not contain duplicates"))
    if !isempty(exps)
        @boundscheck @inbounds exps[begin] ≥ firstindex(e) || throw(BoundsError(e, exps[begin]))
        e isa AbstractExponentsDegreeBounded &&
            (@boundscheck @inbounds exps[end] ≤ lastindex(e) || throw(BoundsError(e, exps[end])))
    end
    return SimpleMonomialVector{Nr,Nc}(unsafe, e, exps)
end

SimpleMonomialVector{Nr,Nc}(exponents::AbstractMatrix{<:Integer}, args...) where {Nr,Nc} =
    SimpleMonomialVector{Nr,Nc,UInt}(exponents, args...)
SimpleMonomialVector{Nr,Nc,I}(exponents::AbstractMatrix{<:Integer}, args...) where {Nr,Nc,I<:Integer} =
    SimpleMonomialVector{Nr,Nc}(ExponentsAll{Nr+2Nc,I}(), exponents, args...)

"""
    SimpleMonomialVector[{I}](mv::AbstractVector{<:AbstractMonomialLike}, along...;
        vars=variables(mv))

Creates a `SimpleMonomialVector` from a generic monomial vector that supports `MultivariatePolynomials`'s interface.
The monomials will internally be represented by the type `I` (`UInt` by default).
The keyword argument `vars` must contain all real-valued and original complex-valued (so not the conjugates) variables that
occur in the monomial vector. However, the order of this iterable (which must have a length) controls how the MP variables are
mapped to [`SimpleVariable`](@ref)s.
The input must not contain duplicates. It will be sorted; if `along` are present, those vectors will be put in the same order
as the inputs.
"""
function SimpleMonomialVector{I}(mv::AbstractVector{<:AbstractMonomialLike}, along::AbstractVector...;
    vars=unique!((x -> isconj(x) ? conj(x) : x).(variables(mv)))) where {I<:Integer}
    any(isconj, vars) && throw(ArgumentError("The variables must not contain conjuates"))
    allunique(vars) || throw(ArgumentError("Variables must not contain duplicates"))

    vars_real = count(isreal, vars)
    vars_complex = count(∘(!, isreal), vars)
    n = length(mv)

    ea = ExponentsAll{vars_real+2vars_complex,I}()
    exps = Vector{I}(undef, n)
    exponents_tmp = Vector{Int}(undef, vars_real + 2vars_complex)
    @inbounds for (j, m) in enumerate(mv)
        i_real, i_complex = 1, 1
        d = 0
        for v in vars
            if isreal(v)
                d += (exponents_tmp[i_real] = degree(m, v))
                i_real += 1
            else
                d += (exponents_tmp[vars_real+2i_complex-1] = degree(m, v))
                d += (exponents_tmp[vars_real+2i_complex] = degree(m, conj(v)))
                i_complex += 1
            end
        end
        exps[j] = exponents_to_index(ea, exponents_tmp, d)
    end
    sort_along!(exps, along...)
    _sortedallunique(exps) || throw(ArgumentError("Monomial vector must not contain duplicates"))
    return SimpleMonomialVector{vars_real,vars_complex}(unsafe, ea, exps)
end

SimpleMonomialVector(mv::AbstractVector{<:AbstractMonomialLike}, args...; kwargs...) =
    SimpleMonomialVector{UInt}(mv, args...; kwargs...)

const SimpleMonomialVectorComplete{Nr,Nc,I<:Integer,E<:AbstractExponentsDegreeBounded,T<:SimpleMonomial{Nr,Nc,I}} =
    SimpleMonomialVector{Nr,Nc,I,E,T}
const SimpleMonomialVectorSubset{Nr,Nc,I<:Integer,E<:AbstractExponents,T<:SimpleMonomial{Nr,Nc,I}} =
    SimpleMonomialVector{Nr,Nc,I,<:Tuple{E,AbstractVector{I}},T}

function Base.getproperty(x::SimpleMonomialVectorComplete{<:Any,<:Any,I}, f::Symbol) where {I<:Integer}
    if f === :e
        return getfield(x, :data)
    elseif f === :indices
        return one(I):unsafe_cast(I, length(unsafe, getfield(x, :data)))
    else
        return getfield(x, f)
    end
end
function Base.getproperty(x::SimpleMonomialVectorSubset, f::Symbol)
    if f === :e
        return getfield(x, :data)[1]
    elseif f === :indices
        return getfield(x, :data)[2]
    else
        return getfield(x, f) # just to trigger the exception
    end
end

Base.IndexStyle(::Type{<:SimpleMonomialVector}) = IndexLinear()
Base.size(x::SimpleMonomialVector) = (unsafe_cast(Int, length(x.indices)),)
Base.firstindex(::SimpleMonomialVector) = 1
Base.lastindex(x::SimpleMonomialVector) = length(x)
Base.isempty(::SimpleMonomialVectorComplete) = false # AbstractExponents must not be empty
Base.copy(x::SimpleMonomialVectorComplete) = x
Base.copy(x::SimpleMonomialVectorSubset{Nr,Nc}) where {Nr,Nc} = SimpleMonomialVector{Nr,Nc}(unsafe, x.e, copy(x.indices))

@inline function Base.getindex(x::SimpleMonomialVector{Nr,Nc}, i::Integer) where {Nr,Nc}
    @boundscheck checkbounds(x, i)
    @inbounds return SimpleMonomial{Nr,Nc}(unsafe, x.e, x.indices[i])
end
@inline function Base.getindex(x::SimpleMonomialVector{Nr,Nc}, i) where {Nr,Nc}
    @boundscheck checkbounds(x, i)
    @inbounds return SimpleMonomialVector{Nr,Nc}(unsafe, x.e, x.indices[i])
end
@inline function Base.view(x::SimpleMonomialVector{Nr,Nc}, i) where {Nr,Nc}
    @boundscheck checkbounds(x, i)
    @inbounds return SimpleMonomialVector{Nr,Nc}(unsafe, x.e, view(x.indices, i))
end
"""
    keepat!!(x::SimpleMonomialVector, i::AbstractVector{Bool})

Keeps the `j`ᵗʰ monomial in `x` only if `i[j]` is `true`. This will mutate `x` if possible (i.e., if it was already indexed by
a vector before), but it might also create a new vector if required (i.e., if a whole range of exponents was covered). Always
use the return value, never rely on `x`.
This function is not exported.
"""
Base.@propagate_inbounds keepat!!(x::SimpleMonomialVector{Nr,Nc}, i::AbstractVector{Bool}) where {Nr,Nc} =
    SimpleMonomialVector{Nr,Nc}(unsafe, x.e, keepat!!(x.indices, i))
Base.@propagate_inbounds keepat!!(x::Vector, i) = keepat!(x, i)
@inline function keepat!!(r::AbstractRange{T}, i::AbstractVector{Bool}) where {T}
    @boundscheck checkbounds(r, length(i))
    len = count(i)
    result = Vector{T}(undef, len)
    j = 1
    for (rᵢ, take) in zip(r, i)
        if take
            @inbounds result[j] = rᵢ
            j += 1
        end
    end
    return result
end

struct SimpleMonomialVectorIterator{Indexed,V,MV,Iterate}
    v::V
    mv::MV

    function SimpleMonomialVectorIterator{Indexed}(v::V, mv::MV) where {Indexed,V<:AbstractVector{Int},Nr,Nc,MV<:SimpleMonomialVector{Nr,Nc}}
        Indexed isa Bool || throw(MethodError(SimpleMonomialVectorIterator{Indexed}, (v, mv)))
        length(v) == Nr + 2Nc || throw(DimensionMismatch("Vector length does not match variable count"))
        # iteration is always better if we have a complete covering of the exponents or at least an assuredly continuous one
        new{Indexed,V,MV,
            MV isa SimpleMonomialVector{Nr,Nc,<:Integer,<:AbstractExponents} ||
            MV isa SimpleMonomialVector{Nr,Nc,<:Integer,<:Tuple{AbstractExponents,AbstractUnitRange}} ||
            isempty(mv) || 2length(mv) ≥ @inbounds(mv.indices[end] - mv.indices[begin] +1)
           }(v, mv)
    end

    function SimpleMonomialVectorIterator{Indexed}(mv::MV) where {Indexed,Nr,Nc,MV<:SimpleMonomialVector{Nr,Nc}}
        Indexed isa Bool || throw(MethodError(SimpleMonomialVectorIterator{Indexed}, (mv,)))
        new{Indexed,Nothing,MV,
            MV isa SimpleMonomialVector{Nr,Nc,<:Integer,<:AbstractExponents} ||
            MV isa SimpleMonomialVector{Nr,Nc,<:Integer,<:Tuple{AbstractExponents,AbstractUnitRange}} ||
            isempty(mv) || 2length(mv) ≥ @inbounds(mv.indices[end] - mv.indices[begin] +1)
           }(nothing, mv)
    end
end

Base.IteratorSize(::Type{<:SimpleMonomialVectorIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:SimpleMonomialVectorIterator}) = Base.HasEltype()
Base.length(mi::SimpleMonomialVectorIterator) = length(mi.mv)
Base.eltype(::Type{<:SimpleMonomialVectorIterator{<:Any,Nothing}}) = Vector{Int}
Base.eltype(::Type{<:SimpleMonomialVectorIterator{<:Any,V}}) where {V} = V

function Base.iterate(mi::SimpleMonomialVectorIterator{Indexed,V,<:SimpleMonomialVector{<:Any,<:Any,I,E},Iterate}) where {Indexed,V,I<:Integer,E,Iterate}
    isempty(mi.mv) && return nothing
    if V <: AbstractVector{Int}
        @inbounds v = copyto!(mi.v, exponents(mi.mv[begin]))
    else
        @assert(V === Nothing)
        @inbounds v = collect(exponents(mi.mv[begin]))
    end
    @inbounds if E <: AbstractExponents
        return Indexed ? ((one(I), v), (one(I), v)) : (v, v)
    elseif E <: Tuple{AbstractExponents,AbstractUnitRange}
        return Indexed ? ((mi.mv.indices[begin], v), (mi.mv.indices[begin], v, length(mi.mv) -1)) : (v, (v, length(mi.mv) -1))
    elseif Iterate
        return (Indexed ? (mi.mv.indices[begin], v) : v), (v, mi.mv.indices[begin], 2, length(mi.mv) -1)
    else
        return Indexed ? ((mi.mv.indices[begin], v), (v, 2, length(mi.mv) -1)) : (v, (v, 2, length(mi.mv) -1))
    end
end
Base.iterate(mi::SimpleMonomialVectorIterator{false}, state::AbstractVector{Int}) =
    @inbounds iterate!(unsafe, state, mi.mv.e) ? (state, state) : nothing
Base.iterate(mi::SimpleMonomialVectorIterator{true}, (index, state)::Tuple{Integer,AbstractVector{Int}}) =
    @inbounds iterate!(unsafe, state, mi.mv.e) ? ((index + one(index), state), (index + one(index), state)) : nothing
function Base.iterate(mi::SimpleMonomialVectorIterator{false}, (state, remaining)::Tuple{AbstractVector{Int},Int})
    iszero(remaining) && return nothing
    success = iterate!(unsafe, state, mi.mv.e)
    @assert(success)
    return state, (state, remaining -1)
end
function Base.iterate(mi::SimpleMonomialVectorIterator{true}, (index, state, remaining)::Tuple{Integer,AbstractVector{Int},Int})
    iszero(remaining) && return nothing
    success = iterate!(unsafe, state, mi.mv.e)
    @assert(success)
    return (index + one(index), state), (index + one(index), state, remaining -1)
end
function Base.iterate(mi::SimpleMonomialVectorIterator{Indexed}, (state, previndex, internal, remaining)::Tuple{AbstractVector{Int},Integer,Int,Int}) where {Indexed}
    iszero(remaining) && return nothing
    index = mi.mv.indices[internal]
    @inbounds for _ in previndex+one(previndex):index
        previndex += 1
        success = iterate!(unsafe, state, mi.mv.e)
        @assert(success)
    end
    return (Indexed ? (index, state) : state), (state, index, internal +1, remaining -1)
end
function Base.iterate(mi::SimpleMonomialVectorIterator{Indexed}, (state, index, remaining)::Tuple{AbstractVector{Int},Int,Int}) where {Indexed}
    iszero(remaining) && return nothing
    result = copyto!(state, exponents(mi.mv[index]))
    @inbounds return (Indexed ? (mi.mv.indices[index], result) : result), (state, index +1, remaining -1)
end
Base.parent(mi::SimpleMonomialVectorIterator) = mi.mv

"""
    veciter(mv::SimpleMonomialVector[, v::AbstractVector{Int}], indexed::Bool=false)

Creates an iterator over all exponents present in `mv` (see
[`veciter`](@ref veciter(::AbstractExponents, ::AbstractVector{Int})) for `AbstractExponents`). By setting `indexed` to `true`,
this iterator will instead give a tuple similar to `enumerate`, where the first index corresponds to the index of the monomial
in the exponent set (so it does not necessarily start at `1` or have unit step).
For type stability, `indexed` may instead be `Val(false)` or `Val(true)`.
"""
veciter(mv::SimpleMonomialVector, v::AbstractVector{Int}, ::Val{indexed}=Val(false)) where {indexed} =
    SimpleMonomialVectorIterator{indexed}(v, mv)
veciter(mv::SimpleMonomialVector, v::AbstractVector{Int}, indexed::Bool) = veciter(mv, v, Val(indexed))
veciter(mv::SimpleMonomialVector, ::Val{indexed}=Val(false)) where {indexed} = SimpleMonomialVectorIterator{indexed}(mv)
veciter(mv::SimpleMonomialVector, indexed::Bool) = veciter(mv, Val(indexed))

MultivariatePolynomials.mindegree(x::SimpleMonomialVectorComplete) = x.e.mindeg
MultivariatePolynomials.mindegree(x::SimpleMonomialVectorSubset) =
    isempty(x.indices) ? 0 : degree_from_index(unsafe, x.e, @inbounds x.indices[begin])
MultivariatePolynomials.maxdegree(x::SimpleMonomialVectorComplete) = x.e.maxdeg
MultivariatePolynomials.maxdegree(x::SimpleMonomialVectorSubset) =
    isempty(x.indices) ? 0 : degree_from_index(unsafe, x.e, @inbounds x.indices[end])
MultivariatePolynomials.extdegree_complex(x::SimpleMonomialVector{<:Any,0}) = extdegree(x)
MultivariatePolynomials.mindegree_complex(x::SimpleMonomialVector{<:Any,0}) = mindegree(x)
MultivariatePolynomials.maxdegree_complex(x::SimpleMonomialVector{<:Any,0}) = maxdegree(x)
MultivariatePolynomials.exthalfdegree(x::SimpleMonomialVector{<:Any,0}) = div.(extdegree(x), 2, RoundUp)
MultivariatePolynomials.minhalfdegree(x::SimpleMonomialVector{<:Any,0}) = div(mindegree(x), 2, RoundUp)
MultivariatePolynomials.maxhalfdegree(x::SimpleMonomialVector{<:Any,0}) = div(maxdegree(x), 2, RoundUp)
# TODO: is there any more efficient way to obtain the [ext,min,max][half]degree[_complex] information for complex-valued
# vectors than the default iteration over exponents. This should be possible at least partly, though it might also depend on
# whether there is a multideg restriction present or not. We might be able to infer the maximal degree for a range of indices
# and a range exponents, at least if these exponents are to the right.

Base.isreal(::SimpleMonomialVector{<:Any,0}) = true
Base.conj(v::SimpleMonomialVector{<:Any,0}, along...) = v
Base.conj(v::SimpleMonomialVector, along...) = _conj(v, along...)
_conj(v::SimpleMonomialVectorComplete{<:Any,<:Any,<:Integer,<:ExponentsDegree}, along...) = v
function _conj(v::SimpleMonomialVectorComplete{Nr,Nc,I,<:ExponentsMultideg}) where {Nr,Nc,I<:Integer}
    if @view(v.e.minmultideg[Nr+1:2:end]) == @view(v.e.minmultideg[Nr+2:2:end])
        minmultideg = v.e.minmultideg
    else
        minmultideg = similar(v.e.minmultideg)
        @inbounds copyto!(minmultideg, 1, v.e.minmultideg, 1, Nr)
        @inbounds for i in Nr+1:2:Nr+2Nc
            minmultideg[i] = v.e.minmultideg[i+1]
            minmultideg[i+1] = v.e.minmultideg[i]
        end
    end
    if @view(v.e.maxmultideg[Nr+1:2:end]) == @view(v.e.maxmultideg[Nr+2:2:end])
        maxmultideg = v.e.maxmultideg
    else
        maxmultideg = similar(v.e.maxmultideg)
        @inbounds copyto!(maxmultideg, 1, v.e.maxmultideg, 1, Nr)
        @inbounds for i in Nr+1:2:Nr+2Nc
            maxmultideg[i] = v.e.maxmultideg[i+1]
            maxmultideg[i+1] = v.e.maxmultideg[i]
        end
    end
    minmultideg === v.e.minmultideg && maxmultideg === v.e.maxmultideg && return v
    return SimpleMonomialVector{Nr,Nc}(ExponentsMultideg{Nr+2Nc,I}(v.e.mindeg, v.e.maxdeg, minmultideg, maxmultideg))
end
function _conj(v::SimpleMonomialVectorComplete{Nr,Nc,I,<:ExponentsMultideg}, along₁, alongᵣ...) where {Nr,Nc,I<:Integer}
    # just to get along sorted correctly, we have to do all this extra work
    cmv = _conj(SimpleMonomialVector{Nr,Nc}(unsafe, v.e, collect(firstindex(v.e):lastindex(v.e)), along₁, alongᵣ...))
    # but for the return value, we really don't need any indexing
    return SimpleMonomialVector{Nr,Nc}(ExponentsMultideg{Nr+2Nc,I}(cmv.e.mindeg, cmv.e.maxdeg, cmv.e.minmultideg,
        cmv.e.maxmultideg))
end
_conj(v::SimpleMonomialVectorSubset{Nr,Nc,I}, along...) where {Nr,Nc,I<:Integer} =
    SimpleMonomialVector{Nr,Nc}(v.e, sort_along!(map(x -> conj(x).index, v), along...)[1])
function _conj(v::SimpleMonomialVectorSubset{Nr,Nc,I,ExponentsMultideg}, along...) where {Nr,Nc,I<:Integer}
    if @view(v.e.minmultideg[Nr+1:2:end]) == @view(v.e.minmultideg[Nr+2:2:end])
        minmultideg = v.e.minmultideg
    else
        minmultideg = similar(v.e.minmultideg)
        @inbounds copyto!(minmultideg, 1, v.e.minmultideg, 1, Nr)
        @inbounds for i in Nr+1:2:Nr+2Nc
            minmultideg[i] = v.e.minmultideg[i+1]
            minmultideg[i+1] = v.e.minmultideg[i]
        end
    end
    if @view(v.e.maxmultideg[Nr+1:2:end]) == @view(v.e.maxmultideg[Nr+2:2:end])
        maxmultideg = v.e.maxmultideg
    else
        maxmultideg = similar(v.e.maxmultideg)
        @inbounds copyto!(maxmultideg, 1, v.e.maxmultideg, 1, Nr)
        @inbounds for i in Nr+1:2:Nr+2Nc
            maxmultideg[i] = v.e.maxmultideg[i+1]
            maxmultideg[i+1] = v.e.maxmultideg[i]
        end
    end
    enew = minmultideg === v.e.minmultideg && maxmultideg === v.e.maxmultideg ? v.e :
        ExponentsMultideg{Nr+2Nc,I}(v.e.mindeg, v.e.maxdeg, minmultideg, maxmultideg)
    return SimpleMonomialVector{Nr,Nc}(
        enew,
        sort_along!(map(let enew=enew; x -> exponents_to_index(enew, exponents(SimpleConjMonomial(x)), degree(x)) end, v),
            along...)[1]
    )
end

Base.in(x::SimpleMonomial{Nr,Nc}, v::SimpleMonomialVectorComplete{Nr,Nc}) where {Nr,Nc} =
    x.e == v.e || !iszero(convert_index(v.e, x.e, x.index, degree(x)))
Base.in(x::SimpleMonomial{Nr,Nc}, v::SimpleMonomialVectorSubset{Nr,Nc}) where {Nr,Nc} =
    insorted(monomial_index(v.e, x), v.indices)
Base.hasfastin(::Type{<:SimpleMonomialVector}) = true
Base.searchsortedlast(v::SimpleMonomialVectorComplete{Nr,Nc}, x::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc} =
    monomial_index(v.e, x)
Base.searchsortedlast(v::SimpleMonomialVectorSubset{Nr,Nc}, x::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc} =
    searchsortedlast(v.indices, monomial_index(v.e, x))
function Base.searchsortedfirst(v::SimpleMonomialVectorComplete{Nr,Nc}, x::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc}
    i = monomial_index(v.e, x)
    return iszero(i) ? lastindex(v) +1 : i
end
function Base.searchsortedfirst(v::SimpleMonomialVectorSubset{Nr,Nc}, x::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc}
    i = monomial_index(v.e, x)
    iszero(i) && return lastindex(v) +1
    return searchsortedfirst(v.indices, i)
end
function Base.searchsorted(v::SimpleMonomialVectorComplete{Nr,Nc}, x::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc}
    i = monomial_index(v.e, x)
    iszero(i) && return 1:0
    return i:i
end
function Base.searchsorted(v::SimpleMonomialVectorSubset{Nr,Nc}, x::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc}
    i = monomial_index(v.e, x)
    iszero(i) && return 1:0
    idx = searchsortedlast(v.indices, i)
    iszero(idx) && return 1:0
    return idx:idx
end
function _issubset_generic(a::SimpleMonomialVector{Nr,Nc}, b::SimpleMonomialVector{Nr,Nc}) where {Nr,Nc}
    isempty(a) && return true
    length(a) > length(b) && return false
    @inbounds if a.e == b.e
        ai = a.indices
        bi = let bi=b.indices
            firsti = searchsortedfirst(bi, ai[1])
            (firsti > lastindex(bi) || bi[firsti] != ai[1]) && return false
            lasti = searchsortedlast(bi, ai[end])
            (lasti < firstindex(bi) || bi[lasti] != ai[end]) && return false
            length(ai) > lasti - firsti +1 && return false
            length(ai) ≤ 2 && return true
            @view(bi[firsti:lasti])
        end
        la, lb = lastindex(ai), lastindex(bi)
        ia = 2
        ib = searchsortedfirst(bi, ai[ia])
        ib > lb && return false
        while true
            while ai[ia] == bi[ib]
                ia += 1
                ia ≥ la && return true
                ib += 1 # we cannot exceed bi, as the last entry is present
            end
            if ai[ia] > bi[ib]
                ib = ib + searchsortedfirst(@view(bi[ib+1:end]), ai[ia])
                ib > lb && return false
            else
                return false
            end
        end
    else
        ae, be = a.e, b.e
        ai = a.indices
        bi = let bi=b.indices
            va = convert_index(be, ae, ai[1])
            firsti = searchsortedfirst(bi, va)
            (firsti > lastindex(bi) || bi[firsti] != va) && return false
            va = convert_index(be, ae, ai[end])
            lasti = searchsortedlast(bi, va)
            (lasti < firstindex(bi) || bi[lasti] != va) && return false
            length(ai) > lasti - firsti +1 && return false
            length(ai) ≤ 2 && return true
            @view(bi[firsti:lasti])
        end
        la, lb = lastindex(ai), lastindex(bi)
        ia = 2
        va = convert_index(unsafe, be, ae, ai[ia])
        ib = searchsortedfirst(bi, va)
        ib > lb && return false
        while true
            while va == bi[ib]
                ia += 1
                ia ≥ la && return true
                va = convert_index(unsafe, be, ae, ai[ia])
                ib += 1
            end
            if va > bi[ib]
                ib = ib + searchsortedfirst(@view(bi[ib+1:end]), va)
                ib > lb && return false
            else
                return false
            end
        end
    end
end
function Base.issubset(a::SimpleMonomialVector{Nr,Nc},
    b::SimpleMonomialVector{Nr,Nc,<:Integer,<:Union{<:AbstractExponents,<:Tuple{AbstractExponents,AbstractUnitRange}}}) where {Nr,Nc}
    isempty(a) && return true
    if a.e isa ExponentsMultideg
        if b.e isa ExponentsMultideg
            all(splat(≥), zip(a.e.minmultideg, b.e.minmultideg)) && all(splat(≤), zip(a.e.maxmultideg, b.e.maxmultideg)) &&
                @goto simple
        else
            @goto simple
        end
    elseif b.e isa ExponentsAll || b.e isa ExponentsDegree
        @goto simple
    else
        @assert(b.e isa ExponentsMultideg)
        all(isequal(0), b.e.minmultideg) && all(Base.Fix2(≥, degree(@inbounds a[end])), b.e.maxmultideg) && @goto simple
    end
    return _issubset_generic(a, b) # invoke will not be able to infer the correct method at compile time
    @label simple
    @inbounds return a[begin] ∈ b && a[end] ∈ b
end
Base.issubset(a::SimpleMonomialVectorComplete{Nr,Nc,<:Integer,<:ExponentsDegree},
    b::SimpleMonomialVectorComplete{Nr,Nc,<:Integer,<:ExponentsDegree}) where {Nr,Nc} =
    return a.e.mindeg ≥ b.e.mindeg && a.e.maxdeg ≤ b.e.maxdeg
Base.issubset(a::SimpleMonomialVector{Nr,Nc}, b::SimpleMonomialVector{Nr,Nc}) where {Nr,Nc} = _issubset_generic(a, b)
function Base.:(==)(a::SimpleMonomialVector{Nr,Nc,<:Integer,<:Union{<:AbstractExponents,<:Tuple{AbstractExponents,AbstractUnitRange}}},
    b::SimpleMonomialVector{Nr,Nc,<:Integer,<:Union{<:AbstractExponents,<:Tuple{AbstractExponents,AbstractUnitRange}}}) where {Nr,Nc}
    length(a) == length(b) || return false
    isempty(a) && return true
    @inbounds if a.e == b.e
        return a[begin] == b[begin] && a[end] == b[end]
    else
        return all(splat(isequal), zip(a, b))
    end
end
Base.:(==)(a::SimpleMonomialVectorComplete{Nr,Nc,<:Integer,E},
    b::SimpleMonomialVectorComplete{Nr,Nc,<:Integer,E}) where {Nr,Nc,E<:AbstractExponents} =
    return a.e == b.e
function Base.:(==)(a::SimpleMonomialVector{Nr,Nc}, b::SimpleMonomialVector{Nr,Nc}) where {Nr,Nc}
    length(a) == length(b) || return false
    isempty(a) && return true
    @inbounds if a.e == b.e
        return all(splat(isequal), zip(a.indices, b.indices))
    else
        return all(splat(isequal), zip(a, b))
    end
end
Base.sort!(x::SimpleMonomialVector) = x # we are already sorted; and no keyword argument that changes the order is allowed

MultivariatePolynomials.variables(::XorTX{<:AbstractVector{<:SimpleMonomial{Nr,Nc}}}) where {Nr,Nc} = SimpleVariables{Nr,Nc}()
MultivariatePolynomials.nvariables(::XorTX{AbstractVector{<:SimpleMonomial{Nr,Nc}}}) where {Nr,Nc} = Nr + 2Nc

struct SimpleMonomialVectorEffectiveVariables{MV<:SimpleMonomialVector}
    mv::MV
end

Base.getproperty(smvev::SimpleMonomialVectorEffectiveVariables, f::Symbol) = getproperty(getfield(smvev, :mv), f)

Base.IteratorSize(::Type{<:SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVectorComplete}}) = Base.HasLength()
Base.IteratorSize(::Type{<:SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVectorSubset}}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:SimpleMonomialVectorEffectiveVariables}) = Base.HasEltype()
Base.eltype(::Type{<:SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{Nr,Nc}}}) where {Nr,Nc} =
    SimpleVariable{Nr,Nc,smallest_unsigned(Nr + 2Nc)}

Base.length(::SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{Nr,Nc,<:Integer,<:ExponentsAll}}) where {Nr,Nc} =
    Nr + 2Nc
Base.length(smvev::SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{Nr,Nc,<:Integer,<:ExponentsDegree}}) where {Nr,Nc} =
    iszero(smvev.e.maxdeg) ? zero(Nr + 2Nc) : Nr + 2Nc
function Base.length(smvev::SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{<:Any,<:Any,<:Integer,<:ExponentsMultideg}})
    i = 0
    e = smvev.e
    for (min, max) in zip(e.minmultideg, e.maxmultideg)
        if !iszero(max) && e.Σminmultideg - min < e.maxdeg
            i += 1
        end
    end
    return i
end

function Base.iterate(::SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{Nr,Nc,<:Integer,<:ExponentsAll}}) where {Nr,Nc}
    iszero(Nr + Nc) && return nothing # is Nr = Nc = 0 allowed at all?
    one_ = one(smallest_unsigned(Nr + 2Nc))
    return SimpleVariable{Nr,Nc}(one_), one_ + one_
end
function Base.iterate(smvev::SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{Nr,Nc,<:Integer,<:ExponentsDegree}}) where {Nr,Nc}
    (iszero(Nr + Nc) || iszero(smvev.e.maxdeg)) && return nothing
    one_ = one(smallest_unsigned(Nr + 2Nc))
    return SimpleVariable{Nr,Nc}(one_), one_ + one_
end
function Base.iterate(::SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{Nr,Nc,<:Integer,<:Union{<:ExponentsAll,<:ExponentsDegree}}},
    var) where {Nr,Nc}
    var > Nr + 2Nc && return nothing
    return SimpleVariable{Nr,Nc}(var), var + one(var)
end
function Base.iterate(smvev::SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{Nr,Nc,<:Integer,<:ExponentsMultideg}},
    var=one(smallest_unsigned(Nr + 2Nc))) where {Nr,Nc}
    e = smvev.e
    @inbounds while var ≤ Nr + 2Nc && (iszero(e.maxmultideg[var]) || e.Σminmultideg - e.minmultideg[var] ≥ e.maxdeg)
        var += one(var)
    end
    var > Nr + 2Nc && return nothing
    return SimpleVariable{Nr,Nc}(var), var + one(var)
end
function Base.iterate(smvev::SimpleMonomialVectorEffectiveVariables{<:SimpleMonomialVector{Nr,Nc,I}},
    var=one(smallest_unsigned(Nr + 2Nc))) where {Nr,Nc,I<:Integer}
    (var > Nr + 2Nc || iszero(var)) && return nothing # zero might happen if var + one(var) overflowed previously
    e = smvev.e
    indices = smvev.indices
    isempty(indices) && return nothing
    @inbounds if e isa ExponentsMultideg
        while iszero(e.maxmultideg[var]) || e.Σminmultideg - e.minmultideg[var] ≥ e.maxdeg
            var += one(var)
            var > Nr + 2Nc && return nothing
        end
        iszero(e.minmultideg[var]) || return SimpleVariable{Nr,Nc}(var), var + one(var)
    end
    counts = index_counts(unsafe, e)
    # This implementation is the zero-allocation variant: go through the indices, check whether the variable is present; if
    # not, skip over as many as possible. As var increases, the size of the skips decreases, so this becomes more and more
    # inefficient. An alternative approach works differently by checking whether indices are present. This requires to create
    # an ExponentsMultideg (with cache) that needs the given variable in its minmultideg.
    # - If indices is an AbstractUnitStep, do a binary search through this multideg to find the first index that is present in
    #   the range.
    # - If indices is anything else, get the first item and (if it is not present), do an exponential search in order to find
    #   the most distant monomial in the multideg whose distance in e is the same (this is the range that contains monomials
    #   where the variable is set). If the intersection is empty, continue with the next item.
    while true
        i = 1
        @inbounds while i ≤ length(indices)
            index = indices[i]
            exps = exponents_from_index(unsafe, e, index)
            lastsameidx, lastexp = exponents_to_index(e, exps, sum(exps), var)
            @assert(!iszero(lastsameidx) && lastexp ≥ 0)
            # now, exponents_to_index only got the first var entries of the exponents. The last of this elements is stored in
            # exps_iter.last.
            # Even if the exponent is nonzero, this does not mean that it occurs, as we only checked the first elements.
            if !iszero(lastexp)
                # the exponent occurs. But this now requires that there is an index in [index, lastsameindex] that can be found
                # in indices.
                return SimpleVariable{Nr,Nc}(var), var + one(var)
            end
            # the variable is zero. lastsameidx is the largest index that has the current set of exponents up to var. Increment
            # by one and find the next available index.
            index = lastsameidx + one(lastsameidx)
            i += searchsortedfirst(@view(indices[i+1:end]), index) # TODO: an exponential search would be better, we can expect
                                                                   # this to be shifted to the start
        end
        var == Nr + 2Nc && return nothing
        var += one(var)
    end
end

MultivariatePolynomials.effective_variables(x::SimpleMonomialVector) = SimpleMonomialVectorEffectiveVariables(x)

"""
    effective_nvariables(x::Union{<:SimpleMonomialVector{Nr,Nc},
                                  <:AbstractArray{<:SimpleMonomialVector{Nr,Nc}}}...)

Calculates the number of effective variable of its arguments: there are at most `Nr + 2Nc` variables that may occur in any
of the monomial vectors or arrays of monomial vectors in the arguments. This function calculates efficiently the number of
variables that actually occur at least once anywhere in any argument.
"""
effective_nvariables(x::SimpleMonomialVectorComplete) = length(effective_variables(x))
# TODO (maybe): get rid off this generated function, instead define an iterable that iterates through all effective variables
# of products/array elements. Maybe also convert every index to ExponentsAll unless they all have the same exponent set, so
# that we don't have to do expensive comparisons over an over...
@generated function effective_nvariables(x::Union{<:SimpleMonomialVector{Nr,Nc},<:AbstractArray{<:SimpleMonomialVector{Nr,Nc}}}...) where {Nr,Nc}
    items = length(x)
    result = Expr(:block, :(result = 0), :(minitem = nothing))
    varloop = Expr(:block, :(result += 1))
    varloop_end = Expr(:block, :(minitem = nothing))
    for (k, xₖ) in enumerate(x)
        xef = Symbol(:xef, k)
        iter = Symbol(:iter, k)
        if xₖ <: SimpleMonomialVector
            push!(result.args, quote
                $xef = effective_variables(x[$k])
                $iter = iterate($xef)
            end)
            push!(varloop.args, quote
                if !isnothing($iter) && $iter[1] == minitem
                    $iter = iterate($xef, $iter[2])
                end
            end)
            push!(varloop_end.args, quote
                if !isnothing($iter) && (isnothing(minitem) || $iter[1] < minitem)
                    minitem = $iter[1]
                end
            end)
        else
            push!(result.args, quote
                $xef = effective_variables.(x[$k])
                $iter = similar($xef, Base.promote_op(iterate, eltype($xef))) # necessary so that broadcasting doesn't
                $iter .= iterate.($xef)                                       # remove the Nothing capability
            end)
            push!(varloop.args, quote
                for (i, item) in enumerate($xef)
                    if !isnothing($iter[i]) && $iter[i][1] == minitem
                        $iter[i] = iterate(item, $iter[i][2])
                    end
                end
            end)
            push!(varloop_end.args, quote
                for iter in $iter
                    if !isnothing(iter) && (isnothing(minitem) || iter[1] < minitem)
                        minitem = iter[1]
                    end
                end
            end)
        end
    end
    push!(result.args, quote
        while true
            $varloop_end
            isnothing(minitem) && return result
            $varloop
        end
    end)
    return :(@inbounds($result))
end

"""
    monomials(Nr, Nc, degree::AbstractUnitRange{<:Integer};
        minmultideg=nothing, maxmultideg=nothing, filter_exps=nothing,
        filter_mons=nothing, I=UInt)

Returns a [`SimpleMonomialVector`](@ref) with `Nr` real and `Nc` complex variables, total degrees contained in `degree`,
ordered according to `Graded{LexOrder}` and individual variable degrees varying between `minmultideg` and `maxmultideg` (where
real variables come first, then complex variables, then their conjugates).

The monomial vector will take possession of the min/maxmultidegs, do not modify them afterwards.

An additional filter may be employed to drop monomials during the construction. Note that the presence of a filter function
will change to a less efficient internal representation. The filter function can get a vector with the exponents as its
argument (`filter_exps`) or a filter function that gets the corresponding `SimpleMonomial`. The former is more efficient if
every exponent has to be retrieved (but do not alter the argument).

The internal representation will be of the type `I`.

This function can be made type-stable by passing `Nr` and `Nc` as `Val`s.
"""
function MultivariatePolynomials.monomials(::Val{Nr}, ::Val{Nc}, degree::AbstractUnitRange{<:Integer};
    minmultideg::Union{Nothing,<:AbstractVector{<:Integer}}=nothing,
    maxmultideg::Union{Nothing,<:AbstractVector{<:Integer}}=nothing,
    filter_exps=nothing, filter_mons=nothing, I::Type{<:Integer}=UInt) where {Nr,Nc}
    n = Nr + 2Nc
    if isnothing(minmultideg)
        if isnothing(maxmultideg)
            e = ExponentsDegree{Nr+2Nc,I}(first(degree), last(degree))
        else
            e = ExponentsMultideg{Nr+2Nc,I}(first(degree), last(degree), ConstantVector(0, n), maxmultideg)
        end
    elseif isnothing(maxmultideg)
        e = ExponentsMultideg{Nr+2Nc,I}(first(degree), last(degree), minmultideg, ConstantVector(last(degree), n))
    else
        e = ExponentsMultideg{Nr+2Nc,I}(first(degree), last(degree), minmultideg, maxmultideg)
    end
    isnothing(filter_exps) && isnothing(filter_mons) && return SimpleMonomialVector{Nr,Nc}(e)
    index_counts(e, 0) # populate the cache
    indices = FastVec{I}()
    if !isnothing(filter_exps)
        i = one(I)
        for exps in veciter(e)
            filter_exps(exps) &&
                (isnothing(filter_mons) || filter_mons(SimpleMonomial{Nr,Nc}(unsafe, e, i, sum(exps, init=0)))) &&
                push!(indices, i)
            i += one(I)
        end
    else
        for i in one(I):I(length(e))
            filter_mons(SimpleMonomial{Nr,Nc}(unsafe, e, i)) && push!(indices, i)
        end
    end
    return SimpleMonomialVector{Nr,Nc}(unsafe, e, finish!(indices))
end

MultivariatePolynomials.monomials(Nr::Integer, Nc::Integer, degree::AbstractUnitRange{<:Integer}; kwargs...) =
    monomials(Val(Nr), Val(Nc), degree; kwargs...)

"""
    intersect(a::SimpleMonomialVector{Nr,Nc}, b::SimpleMonomialVector{Nr,Nc})

Calculates efficiently the intersection of two monomial vectors.
"""
Base.intersect(::SimpleMonomialVector{Nr,Nc}, ::SimpleMonomialVector{Nr,Nc}) where {Nr,Nc}

function Base.intersect(a::MV, b::MV) where {Nr,Nc,I<:Integer,MV<:SimpleMonomialVector{Nr,Nc,I,<:ExponentsDegree}}
    a.e == b.e && return a
    ae, be = a.e, b.e
    mindeg = max(ae.mindeg, be.mindeg)
    maxdeg = min(ae.maxdeg, be.maxdeg)
    return mindeg ≤ maxdeg ? SimpleMonomialVector{Nr,Nc}(ExponentsDegree{Nr+2Nc,I}(mindeg, maxdeg)) :
                             SimpleMonomialVector{Nr,Nc}(unsafe, ae, I[])
end
function Base.intersect(a::SimpleMonomialVector{Nr,Nc,I,<:ExponentsMultideg},
    b::SimpleMonomialVector{Nr,Nc,I,<:ExponentsMultideg}) where {Nr,Nc,I<:Integer}
    # potentially type-unstable
    a.e == b.e && return a
    ae, be = a.e, b.e
    mindeg = max(ae.mindeg, be.mindeg)
    maxdeg = min(ae.maxdeg, be.maxdeg)
    mindeg ≤ maxdeg || return SimpleMonomialVector{Nr,Nc}(unsafe, ae, I[])
    if ae.minmultideg == be.minmultideg
        minmultideg = ae.minmultideg
    else
        minmultideg = elementwise(max, ae.minmultideg, be.minmultideg)
    end
    if ae.maxmultideg == be.maxmultideg
        maxmultideg = ae.maxmultideg
    else
        maxmultideg = elementwise(min, ae.maxmultideg, be.maxmultideg)
    end
    all(splat(≤), zip(minmultideg, maxmultideg)) || return SimpleMonomialVector{Nr,Nc}(unsafe, ae, I[])
    any(∘(!, iszero), minmultideg) || any(<(maxdeg), maxmultideg) ||
        return SimpleMonomialVector{Nr,Nc}(ExponentsDegree{Nr+2Nc,I}(mindeg, maxdeg))
    return SimpleMonomialVector{Nr,Nc}(ExponentsMultideg{Nr+2Nc,I}(mindeg, maxdeg, minmultideg, maxmultideg))
end
function Base.intersect(a::SimpleMonomialVector{Nr,Nc,I,<:ExponentsDegree},
    b::SimpleMonomialVector{Nr,Nc,I,<:ExponentsMultideg}) where {Nr,Nc,I<:Integer}
    # type-unstable (but only two choices)
    ae, be = a.e, b.e
    ae.mindeg == be.mindeg && ae.maxdeg == be.maxdeg && return b
    mindeg = max(ae.mindeg, be.mindeg)
    maxdeg = min(ae.maxdeg, be.maxdeg)
    mindeg ≤ maxdeg || return SimpleMonomialVector{Nr,Nc}(unsafe, be, I[])
    any(∘(!, iszero), be.minmultideg) || any(<(maxdeg), be.maxmultideg) ||
        return SimpleMonomialVector{Nr,Nc}(ExponentsDegree{Nr+2Nc,I}(mindeg, maxdeg))
    return SimpleMonomialVector{Nr,Nc}(ExponentsMultideg{Nr+2Nc,I}(mindeg, maxdeg, be.minmultideg, be.maxmultideg))
end
Base.intersect(a::SimpleMonomialVector{Nr,Nc,I,<:ExponentsMultideg},
    b::SimpleMonomialVector{Nr,Nc,I,<:ExponentsDegree}) where {Nr,Nc,I<:Integer} = intersect(b, a)
function Base.intersect(a::SimpleMonomialVector{Nr,Nc,I,<:ExponentsDegree},
    b::SimpleMonomialVectorSubset{Nr,Nc,I}) where {Nr,Nc,I<:Integer}
    @assert(!isempty(a))
    isempty(b) && return SimpleMonomialVector{Nr,Nc}(unsafe, be, @view(b.indices[begin:end])) # just to have type-stability
    ae, be = a.e, b.e
    @inbounds if ae.mindeg ≤ degree(b[begin])
        if ae.maxdeg ≥ degree(b[end])
            return SimpleMonomialVector{Nr,Nc}(unsafe, be, @view(b.indices[begin:end]))
        else
            lastina = searchsortedlast(b.indices, convert_index(unsafe, be, ae, lastindex(ae), ae.maxdeg))
        end
        firstina = 1
    elseif ae.maxdeg ≥ degree(b[end])
        firstina = searchsortedfirst(b.indices, convert_index(unsafe, be, ae, firstindex(ae), ae.mindeg))
        lastina = length(b.indices)
    else
        firstina = searchsortedfirst(b.indices, convert_index(unsafe, be, ae, firstindex(ae), ae.mindeg))
        lastina = searchsortedlast(b.indices, convert_index(unsafe, be, ae, lastindex(ae), ae.maxdeg))
    end
    return SimpleMonomialVector{Nr,Nc}(unsafe, be, @view(b.indices[firstina:lastina]))
end
Base.intersect(a::SimpleMonomialVectorSubset{Nr,Nc,I},
    b::SimpleMonomialVector{Nr,Nc,I,<:ExponentsDegree}) where {Nr,Nc,I<:Integer} = intersect(b, a)
function Base.intersect(a::Union{<:SimpleMonomialVector{Nr,Nc,I,<:ExponentsMultideg},
                                 <:SimpleMonomialVectorSubset{Nr,Nc,I}}, b::SimpleMonomialVectorSubset{Nr,Nc,I}) where {Nr,Nc,I<:Integer}
    # necessarily type-unstable (but finite Union)
    ae, be = a.e, b.e
    isempty(a) && return a
    isempty(b) && return b
    if ae == be
        return SimpleMonomialVector{Nr,Nc}(unsafe, ae, intersect_sorted(a.indices, b.indices))
    else
        @assert(!(ae isa ExponentsAll && be isa ExponentsAll))
        mindeg = max(ae isa AbstractExponentsDegreeBounded ? ae.mindeg : 0,
                     be isa AbstractExponentsDegreeBounded ? be.mindeg : 0)
        maxdeg = min(ae isa AbstractExponentsDegreeBounded ? ae.maxdeg : typemax(Int),
                     be isa AbstractExponentsDegreeBounded ? be.maxdeg : typemax(Int))
        mindeg ≤ maxdeg || return SimpleMonomialVector{Nr,Nc}(unsafe, ae, I[])
        if ae isa ExponentsMultideg || be isa ExponentsMultideg
            if !(ae isa ExponentsMultideg)
                minmultideg = be.minmultideg
                maxmultideg = be.maxmultideg
            elseif !(be isa ExponentsMultideg)
                minmultideg = ae.minmultideg
                maxmultideg = ae.maxmultideg
            else
                minmultideg = elementwise(max, ae.minmultideg, be.minmultideg)
                maxmultideg = elementwise(min, ae.maxmultideg, be.maxmultideg)
                all(splat(≤), zip(minmultideg, maxmultideg)) || return SimpleMonomialVector{Nr,Nc}(unsafe, ae, I[])
            end
            newe = ExponentsMultideg{Nr+2Nc,I}(mindeg, maxdeg, minmultideg, maxmultideg)
        else
            newe = ExponentsDegree{Nr+2Nc,I}(mindeg, maxdeg)
        end
        rema, remb = length(a), length(b)
        (iszero(rema) || iszero(remb)) && return SimpleMonomialVector{Nr,Nc}(unsafe, newe, I[])
        index_counts(newe, 0) # populate the cache
        indices_a = a.indices
        indices_b = b.indices
        indices = FastVec{I}(buffer=length(newe))
        ia, ib = 1, 1
        # TODO (maybe): there's a potentially better version for vectors with small overlap. We can convert the smaller index
        # to the other exponent set, do a binary (or exponential) search in this set, then adjust the index. This allows to
        # very quickly skip over all nonintersecting ranges. When we found equality, we need to convert to the new index set.
        # However, when the overlap is large, this introduces unnecessary overhead.
        @inbounds begin
            nexta = convert_index(unsafe, newe, ae, indices_a[ia])
            nextb = convert_index(unsafe, newe, be, indices_b[ib])
            while true
                while iszero(nexta) || nexta < nextb
                    iszero(rema -= 1) && @goto done
                    ia += 1
                    nexta = convert_index(unsafe, newe, ae, indices_a[ia])
                end
                while iszero(nextb) || nextb < nexta
                    iszero(remb -= 1) && @goto done
                    ib += 1
                    nextb = convert_index(unsafe, newe, be, indices_b[ib])
                end
                if nexta == nextb
                    unsafe_push!(indices, nexta)
                    iszero(rema -= 1) && @goto done
                    iszero(remb -= 1) && @goto done
                    ia += 1
                    nexta = convert_index(unsafe, newe, ae, indices_a[ia])
                    ib += 1
                    nextb = convert_index(unsafe, newe, be, indices_b[ib])
                end
            end
        end
        @label done
        return SimpleMonomialVector{Nr,Nc}(unsafe, newe, finish!(indices))
    end
end
Base.intersect(a::SimpleMonomialVectorSubset{Nr,Nc,I},
    b::SimpleMonomialVector{Nr,Nc,I,<:ExponentsMultideg}) where {Nr,Nc,I<:Integer} = intersect(b, a)

"""
    merge_monomial_vectors(::Val{Nr}, ::Val{Nc}, e::AbstractExponents, X::AbstractVector)

Returns the vector of monomials contained `X` in increasing order and without any duplicates. The individual elements in
`X` must be sorted iterables with a length and return `SimpleMonomial`s compatible with the number of real `Nr` and complex
variables `Nc`. The output will internally use the exponents `e`.
The result type will be a `SimpleMonomialVector` with `Nr` and `Nc` as given, `I` and the exponents determined by `e`, and
indexed internally with a `Vector{I}`.
"""
function MultivariatePolynomials.merge_monomial_vectors(::Val{Nr}, ::Val{Nc}, e::AbstractExponents{N,I},
    X::AbstractVector) where {Nr,Nc,N,I<:Integer}
    @nospecialize e X
    N == Nr + 2Nc || throw(MethodError(merge_monomial_vectors, (Val(Nr), Val(Nc), e, X)))
    # We must be very careful here. The function must be fast, there's no way to allow for dynamic dispatch within the loop.
    # However, we might get a collection of different types in the the vector (different AbstractExponents, indexed and full).
    # It is certainly not viable to compile a new function for every possible combination of inputs. Therefore, here we first
    # sort them all by type, then pass them to a generated function. But we can bypass this if there is really just one known
    # concrete eltype.
    if isconcretetype(eltype(X))
        return merge_monomial_vectors_impl(Val(Nr), Val(Nc), e, X)
    else
        grouped = Dict{DataType,Vector{<:AbstractVector}}()
        for Xᵢ in X
            T = typeof(Xᵢ)
            if Base.IteratorEltype(T) isa Base.HasEltype
                eltype(T) <: SimpleMonomial{Nr,Nc} || SimpleMonomial{Nr,Nc} <: eltype(T) ||
                    error("An iterator does not have a compatible element type")
            end
            v = get!(let T=T; () -> T[] end, grouped, T)
            push!(v, Xᵢ)
        end
        # We need to sort the types in an arbitrary, but consistent manner. The function is commutative, so no need to generate
        # two different functions just because the order differs. Assuming hash probably never collides on the few couple of
        # types that are possible, we'll use this as a comparison between Type objects.
        return merge_monomial_vectors_impl(Val(Nr), Val(Nc), e,
            sort_along!(hash.(keys(grouped)), collect(values(grouped)))[2]...)
    end
end

_extract_n(::XorTX{<:SimpleMonomial{Nr,Nc}}) where {Nr,Nc} = (Nr, Nc)
_extract_n(_) = missing
_extract_E(::Type{<:SimpleMonomial{<:Any,<:Any,<:Integer,E}}) where {E<:AbstractExponents} = E
_extract_E(_) = Missing
_extract_I(::XorTX{<:AbstractExponents{<:Any,I}}) where {I<:Integer} = I
_extract_I(_) = UInt
_is_full_smv(::XorTX{<:SimpleMonomialVector{<:Any,<:Any,<:Integer,<:AbstractExponents}}) = true
_is_full_smv(_) = false
_to_I(::Type{I}, e::AbstractExponents{<:Any,I}) where {I<:Integer} = e
_to_I(::Type{I}, e::ExponentsAll{N,<:Integer}) where {I<:Integer,N} = ExponentsAll{N,I}()
_to_I(::Type{I}, e::ExponentsDegree{N,<:Integer}) where {I<:Integer,N} = ExponentsDegree{N,I}(e.mindeg, e.maxdeg)
_to_I(::Type{I}, e::ExponentsMultideg{N,<:Integer}) where {I<:Integer,N} =
    ExponentsMultideg{N,I}(e.mindeg, e.maxdeg, e.minmultideg, e.maxmultideg)

@doc """
    merge_monomial_vectors(X::AbstractVector)

Potentially type-unstable variant that automatically determines the output format.
If `X` has a defined eltype with known eltype `<:SimpleMonomial{Nr,Nc,I,E}`, `Nr`, `Nc`, and `I` are determined automatically
in a type-stable manner. If not, they are taken from the eltype of the first iterable in `X` (which is not type stable). If
this is not possible, either, they are taken from the first element in the first nonempty iterable in `X`.

Regarding the automatic determination of `E`, the following rule is applied: it is assumed that if `E` is known in the element
type, then every monomial will have the same instance of `e` as the exponents _per iterable_ (this is always satisfied if the
iterables are `SimpleMonomialVector`s). If there is one exponent that covers all others, this instance will be used. If not,
the largest necessary exponents will be constructed; the result will be indexed unless all can be merged contiguously). Note
that inhomogeneous iterables must implement `last` if the elements are based on `ExponentsAll`.

!!! info
    This method has a very general type signature and may therefore also be called for other implementations of
    `MultivariatePolynomials`. However, this case will be caught and then forwarded to the generic MP implementation.
"""
MultivariatePolynomials.merge_monomial_vectors(X::AbstractVector{<:SimpleMonomialVector})
# ^ Julia assigns docstrings in a faulty way when type parameters are involved. As there is the generic (::Any) function and
# every (unconstrained) type parameter that is present will lead to a pseudo-signature (::Any), we must avoid this.
# (https://github.com/JuliaLang/julia/pull/53824)

function MultivariatePolynomials.merge_monomial_vectors(X::AbstractVector{T}) where {T}
    if Base.IteratorEltype(T) isa Base.HasEltype
        nrnc = _extract_n(eltype(T))
        E = _extract_E(eltype(T))
    else
        nrnc = missing
        E = Missing
    end
    if ismissing(nrnc)
        for Xᵢ in X
            item = iterate(Xᵢ)
            if !isnothing(item)
                if ismissing((nrnc = _extract_n(item[1]);))
                    # oops... Looks like this was not really SimplePolynomials related. Call the generic MP version.
                    return @invoke merge_monomial_vectors(X::Any)
                end
                break
            end
        end
    end
    ismissing(nrnc) && return []
    Nr, Nc = nrnc
    if E <: ExponentsAll
        e = ExponentsAll{Nr+2Nc,_extract_I(E)}()
        no_indexing = _is_full_smv(T)
    else
        e = missing
        no_indexing = true
        copied = false
        @inbounds for Xᵢ in X
            if Xᵢ isa SimpleMonomialVector
                isempty(Xᵢ) && continue
                eᵢ = Xᵢ.e
                no_indexingᵢ = _is_full_smv(Xᵢ)
            else
                item = iterate(Xᵢ)
                isnothing(item) && continue
                eᵢ = item[1].e
                no_indexingᵢ = false
            end
            if ismissing(e)
                if eᵢ isa ExponentsAll
                    e = ExponentsDegree{Nr+2Nc,_extract_I(e)}(degree(Xᵢ[begin]), degree(Xᵢ[end]))
                    if Base.haslength(Xᵢ)
                        no_indexing = length(Xᵢ) == length(e)
                    end
                else
                    e = eᵢ
                    no_indexing = no_indexingᵢ
                end
                continue
            end
            I = promote_type(_extract_I(e), _extract_I(eᵢ))
            if eᵢ isa ExponentsAll
                @assert(!no_indexingᵢ) # ExponentsAll can never determine a finite-size vector fully to do the merge, convert
                                       # it to a degree-bound index
                eᵢ = ExponentsDegree{Nr+2Nc,I}(degree(Xᵢ[begin]), degree(Xᵢ[end]))
                if Base.haslength(Xᵢ)
                    no_indexingᵢ = length(Xᵢ) == length(eᵢ)
                end
            end
            if eᵢ isa ExponentsDegree || e isa ExponentsDegree
                e_deg, e_other, noind_deg, noind_other = eᵢ isa ExponentsDegree ? (eᵢ, e, no_indexingᵢ, no_indexing) :
                                                                                  (e, eᵢ, no_indexing, no_indexingᵢ)
                if e_other isa ExponentsDegree
                    # are the exponents identical? - we can drop indexing if at least one doesn't need any
                    if e_deg.mindeg == e_other.mindeg && e_deg.maxdeg == e_other.maxdeg
                        if _extract_I(e) !== I
                            e = _to_I(I, eᵢ)
                        end
                        no_indexing = noind_deg | noind_other
                        continue
                    end
                    # do the other exponents cover the degree ones fully? - the other ones remain
                    if e_deg.mindeg ≥ e_other.mindeg && e_deg.maxdeg ≤ e_other.maxdeg
                        e, no_indexing = _to_I(I, e_other), noind_other
                        continue
                    end
                else
                    @assert(e_other isa ExponentsMultideg)
                    # do the other exponents cover the degree ones fully? - the other ones remain
                    if e_deg.mindeg ≥ e_other.mindeg && e_deg.maxdeg ≤ e_other.maxdeg && all(iszero, e_other.minmultideg) &&
                        all(≥(e_deg.maxdeg), e_other.maxmultideg)
                        e, no_indexing = _to_I(I, e_other), noind_other
                        continue
                    end
                end
                # do the degree ones cover the other ones fully? - the degree ones remain
                if e_other.mindeg ≥ e_deg.mindeg && e_other.maxdeg ≤ e_deg.maxdeg
                    e, no_indexing = _to_I(I, e_deg), noind_deg
                else
                    if e_other isa ExponentsMultideg && e_other.maxdeg > e_deg.maxdeg && e_other.mindeg ≤ e_deg.maxdeg +1 &&
                        all(iszero, e_other.minmultideg) && all(≥(e_deg.maxdeg), e_other.maxmultideg)
                        if e_other.mindeg == e_deg.mindeg
                            e = _to_I(I, e_other)
                        else
                            e = ExponentsMultideg{Nr+2Nc,I}(e_deg.mindeg, e_other.maxdeg, e_other.minmultideg,
                                    e_other.maxmultideg)
                        end
                        noind = noind_other
                        continue
                    end
                    # are we overlapping or adjacent?
                    no_indexing = noind_deg && noind_other && e_other isa ExponentsDegree &&
                        (e_deg.mindeg ≤ e_other.mindeg ≤ e_deg.maxdeg +1 || e_other.mindeg ≤ e_deg.mindeg ≤ e_other.maxdeg +1)
                    # the exponents don't include each other. Construct covering ones
                    e = ExponentsDegree{Nr+2Nc,I}(min(e_deg.mindeg, e_other.mindeg), max(e_deg.maxdeg, e_other.maxdeg))
                end
            else
                @assert(eᵢ isa ExponentsMultideg && e isa ExponentsMultideg)
                # are the indices identical?
                if eᵢ.mindeg == e.mindeg && eᵢ.maxdeg == e.maxdeg && eᵢ.minmultideg == e.minmultideg &&
                    eᵢ.maxmultideg == e.maxmultideg
                    if _extract_I(e) !== I
                        e = _to_I(I, eᵢ)
                    end
                    no_indexing |= no_indexingᵢ
                    continue
                end
                # do the existing exponents cover the new ones fully?
                if eᵢ.mindeg ≥ e.mindeg && eᵢ.maxdeg ≤ e.maxdeg && all(splat(≥), zip(eᵢ.minmultideg, e.minmultideg)) &&
                    all(splat(≤), zip(eᵢ.maxmultideg, e.maxmultideg))
                    e = _to_I(I, e)
                    continue
                end
                # do the new ones cover the existing ones fully?
                if e.mindeg ≥ eᵢ.mindeg && e.maxdeg ≤ eᵢ.maxdeg && all(splat(≥), zip(e.minmultideg, eᵢ.minmultideg)) &&
                    all(splat(≤), zip(e.maxmultideg, eᵢ.maxmultideg))
                    e, no_indexing = _to_I(I, eᵢ), no_indexingᵢ
                    continue
                else
                    # are we overlapping or adjacent?
                    no_indexing &= no_indexingᵢ && ((e.mindeg ≤ eᵢ.mindeg ≤ e.maxdeg +1 && e.minmultideg == eᵢ.minmultideg &&
                        all(≥(e.maxdeg), eᵢ.maxmultideg)) || (eᵢ.mindeg ≤ e.mindeg ≤ eᵢ.maxdeg +1 &&
                        e.minmultideg == eᵢ.minmultideg && all(≥(eᵢ.maxdeg), e.maxmultideg)))
                    # the exponents don't include each other. Construct covering ones.
                    if !copied
                        e = ExponentsMultideg{Nr+2Nc,I}(min(e.mindeg, eᵢ.mindeg), max(e.maxdeg, eᵢ.maxdeg),
                                elementwise(min, e.minmultideg, eᵢ.minmultideg),
                                elementwise(max, e.maxmultideg, eᵢ.maxmultideg))
                        copied = true
                    else
                        mind = elementwise!(e.minmultideg, min, e.minmultideg, eᵢ.minmultideg)
                        maxd = elementwise!(e.maxmultideg, max, e.maxmultideg, eᵢ.maxmultideg)
                        e = ExponentsMultideg{Nr+2Nc,I}(min(e.mindeg, eᵢ.mindeg), max(e.maxdeg, eᵢ.maxdeg), mind, maxd)
                    end
                end
            end
        end
        if ismissing(e)
            E <: Missing && return []
            return SimpleMonomialVector{Nr,Nc}(ExponentsAll{Nr+2Nc,_extract_I(E)}(), _extract_I(E)[])
        end
    end
    if no_indexing
        return SimpleMonomialVector{Nr,Nc}(e)
    else
        return merge_monomial_vectors(Val(Nr), Val(Nc), e, X)
    end
end

@generated function merge_monomial_vectors_impl(::Val{Nr}, ::Val{Nc}, e::E, Xs::Vector...) where {Nr,Nc,N,I<:Integer,E<:AbstractExponents{N,I}}
    @assert(Nr + 2Nc == N)
    # Here, every vector contained in a single Xs[i] is of the same type and we can specialize
    types = length(Xs)
    result = Expr(:block, :(remaining = 0), :(indices = FastVec{I}()))
    iters = [Symbol(:iters, i) for i in 1:types]
    idxs = [Symbol(:idxs, i) for i in 1:types]
    mins = [Symbol(:min, i) for i in 1:types]
    for i in 1:types
        push!(result.args, quote
            # we cannot just do iterate.(Xs[i]) - if none are Nothing, the Vector won't allow it.
            $(iters[i]) = Vector{$(Base.promote_op(iterate, eltype(Xs[i])))}(undef, length(Xs[$i]))
            $(idxs[i]) = similar($(iters[i]), I)
            $(mins[i]) = 1
            curminidx = typemax(I)
            for (j, Xⱼ) in enumerate(Xs[$i])
                $(iters[i])[j] = iterval = iterate(Xⱼ)
                if isnothing(iterval)
                    $(idxs[i])[j] = typemax(I)
                else
                    $(idxs[i])[j] = mi = monomial_index(e, iterval[1])
                    if mi < curminidx
                        $(mins[i]) = j
                        curminidx = mi
                    end
                end
            end
        end)
    end
    process_min = Expr(:if)
    process_min_cur = process_min
    for i in 1:types
        if !isone(i)
            process_min_cur = let process_min_next=Expr(:elseif)
                push!(process_min_cur.args, process_min_next)
                process_min_next
            end
        end
        process_min_i = Expr(:block)
        # Probably Julia won't be able to figure out that !isnothing($(iters[i])) actually holds true. This would be a good
        # case for the @ensure/@assume... macro proposal to complement @assert (Julia issue #51729). But we don't have it
        # (yet). On the other hand, Cthulhu seems to suggest that while $(iters[i])[...] is indeed always a Union (regardless
        # of whether we wrap it in isnothing or not), accessing an index will automatically give the correct type again, even
        # here.
        push!(process_min_i.args, quote
            lastidx = $(idxs[i])[$(mins[i])]
            push!(indices, lastidx)
            while true
                $(iters[i])[$(mins[i])] = nextiter = iterate(Xs[$i][$(mins[i])], $(iters[i])[$(mins[i])][2])
                if isnothing(nextiter)
                    $(idxs[i])[$(mins[i])] = typemax(I)
                else
                    $(idxs[i])[$(mins[i])] = monomial_index(e, nextiter[1])
                end
                $(mins[i]) = argmin($(idxs[i]))
                $(idxs[i])[$(mins[i])] == lastidx || break
            end
        end)
        push!(process_min_cur.args, :(curmin == $i), process_min_i)
    end
    push!(process_min_cur.args, Expr(:break))
    push!(result.args, quote
        col = 1
        while true
            curmin = 0
            curminidx = typemax(I) -1 # Let's assume we don't ever encounter the two largest values, then we can safely do some
                                      # additional comparison logic
            $((
                :(if $(idxs[i])[$(mins[i])] < curminidx
                    curmin = $i
                    curminidx = $(idxs[i])[$(mins[i])]
                else
                    while $(idxs[i])[$(mins[i])] == curminidx
                        # duplicate, skip it (doesn't matter whether it is the global minimum, if not it will be a duplicate
                        # later)
                        $(iters[i])[$(mins[i])] = nextiter = iterate(Xs[$i][$(mins[i])], $(iters[i])[$(mins[i])][2])
                        if isnothing(nextiter)
                            $(idxs[i])[$(mins[i])] = typemax(I)
                        else
                            $(idxs[i])[$(mins[i])] = monomial_index(e, nextiter[1])
                        end
                        $(mins[i]) = argmin($(idxs[i]))
                    end
                end)
                for i in 1:types
            )...)
            $process_min
            col += 1
        end
    end)
    push!(result.args,
        :(return SimpleMonomialVector{Nr,Nc}(unsafe, e, finish!(indices)))
    )
    return :(@inbounds($result))
end

struct FakeMonomialVector{S<:SimpleMonomialVector,V,M} <: AbstractVector{M}
    data::S
    real_vars::Vector{V}
    complex_vars::Vector{V}

    function FakeMonomialVector(data::S, real_vars::Vector{V}, complex_vars::Vector{V}) where {Nr,Nc,S<:SimpleMonomialVector{Nr,Nc},V<:AbstractVariable}
        (length(real_vars) == Nr && length(complex_vars) == Nc) || error("Invalid monomial vector construction")
        new{S,V,monomial_type(V)}(data, real_vars, complex_vars)
    end
end

Base.length(fmv::FakeMonomialVector) = length(fmv.data)
Base.size(fmv::FakeMonomialVector) = (length(fmv.data),)
function Base.getindex(fmv::FakeMonomialVector{S,V,M} where {V,S}, x) where {M}
    mon = fmv.data[x]
    isconstant(mon) && return constant_monomial(M)
    exps = exponents(mon)
    expit = iterate(exps)
    i = 1
    havemon = false
    while !isnothing(expit)
        i > length(fmv.real_vars) && break
        expᵢ, expitdata = expit
        if !iszero(expᵢ)
            if !havemon
                @inbounds mon = fmv.real_vars[i] ^ expᵢ
                havemon = true
            else
                @inbounds mon *= fmv.real_vars[i] ^ expᵢ
            end
        end
        i += 1
        expit = iterate(exps, expitdata)
    end
    i = 1
    while !isnothing(expit)
        i > length(fmv.complex_vars) && break
        expᵢ, expitdata = expit
        if !iszero(expᵢ)
            if !havemon
                @inbounds mon = fmv.complex_vars[i] ^ expᵢ
                havemon = true
            else
                @inbounds mon *= fmv.complex_vars[i] ^ expᵢ
            end
        end
        expᵢ, expitdata = iterate(exps, expitdata)::Tuple
        if !iszero(expᵢ)
            if !havemon
                @inbounds mon = conj(fmv.complex_vars[i]) ^ expᵢ
                havemon = true
            else
                @inbounds mon *= conj(fmv.complex_vars[i]) ^ expᵢ
            end
        end
        i += 1
        expit = iterate(exps, expitdata)
    end
    @assert(isnothing(expit))
    return mon
end

"""
    change_backend(mv::SimpleMonomialVector, variable::AbstractVector{<:AbstractVariable})

Changes a `SimpleMonomialVector` into a different implementation of `MultivariatePolynomials`, where the variables are taken
from the given vector in the order as they appear (but keeping real and complex variables distinct).

This conversion is not particularly efficient, as it works with generic implementations.
"""
function change_backend(mv::SimpleMonomialVector{Nr,Nc}, variables::AbstractVector{V}) where {Nr,Nc,V<:AbstractVariable}
    real_vars = similar(variables, 0)
    complex_vars = similar(real_vars)
    for v in variables
        if isreal(v)
            push!(real_vars, v)
        elseif isconj(v)
            vo = conj(v)
            vo ∈ complex_vars || push!(complex_vars, vo)
        else
            push!(complex_vars, v)
        end
    end
    (length(real_vars) == Nr && length(complex_vars) == Nc) || throw(ArgumentError("Incompatible variables"))
    return monomial_vector(FakeMonomialVector(mv, real_vars, complex_vars))
end
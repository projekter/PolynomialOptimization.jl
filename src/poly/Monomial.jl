export SimpleMonomial, SimpleConjMonomial, monomial_produt, monomial_index

"""
    SimpleMonomial{Nr,Nc,I<:Integer,E<:AbstractExponents} <: AbstractMonomial

`SimpleMonomial` represents a monomial. In order to be used together in operations, the number of real-valued variables `Nr`
and the number of complex-valued variables `Nc` are fixed in the type. The monomial is identified according to its index (of
type `I`) in an exponent set of type `E`. This should be an unsigned type, but to allow for `BigInt`, no such restriction is
imposed.
"""
struct SimpleMonomial{Nr,Nc,I<:Integer,E<:AbstractExponents} <: AbstractMonomial
    e::E
    index::I
    degree::Int # because it is required very often, we'll precalculate it

    # internal functions, don't use.
    @inline function SimpleMonomial{Nr,Nc}(::Unsafe, e::E, index::I, degree::Int=degree_from_index(unsafe, e, index)) where
        {Nr,Nc,N,I<:Integer,E<:AbstractExponents{N,I}}
        N == Nr + 2Nc || throw(MethodError(SimpleMonomial{Nr,Nc,E}, (unsafe, index, degree)))
        new{Nr,Nc,I,E}(e, index, degree)
    end
end

Base.parent(m::SimpleMonomial) = m
SimpleMonomial(m::SimpleMonomial) = m

"""
    SimpleMonomial{Nr,0[,I]}([e::AbstractExponents,]
        exponents_real::AbstractVector{<:Integer})
    SimpleMonomial{0,Nc[,I]}([e::AbstractExponents,]
        exponents_complex::AbstractVector{<:Integer},
        exponents_conj::AbstractVector{<:Integer})
    SimpleMonomial{Nr,Nc[,I]}([e::AbstractExponents,]
        exponents_real::AbstractVector{<:Integer},
        exponents_complex::AbstractVector{<:Integer},
        exponents_conj::AbstractVector{<:Integer})

Creates a monomial within an exponent set `e`. If `e` is omitted, `ExponentsAll{Nr+2Nc,UInt}` is chosen by default.
Alternatively, all three methods may also be called with the index type `I` as a third type parameter, omitting `e`, which then
chooses `ExponentsAll{Nr+2Nc,I}` by default.
"""
SimpleMonomial{Nr,Nc}(::AbstractExponents, ::AbstractVector{<:Integer}...) where {Nr,Nc}

function SimpleMonomial{Nr,0}(e::AbstractExponents{Nr,I}, exponents_real::AbstractVector{<:Integer}) where {Nr,I<:Integer}
    length(exponents_real) == Nr || throw(ArgumentError("Requested $Nr real variables, but got $(length(exponents_real))"))
    degree = sum(exponents_real, init=0)
    idx = exponents_to_index(e, exponents_real, degree)
    @boundscheck iszero(idx) && throw(BoundsError(e, exponents_real))
    return SimpleMonomial{Nr,0}(unsafe, e, idx, degree)
end

function SimpleMonomial{0,Nc}(e::AbstractExponents{N,I}, exponents_complex::AbstractVector{<:Integer},
    exponents_conj::AbstractVector{<:Integer}) where {Nc,N,I<:Integer}
    N == 2Nc || throw(MethodError(SimpleMonomial{0,Nc}, (e, exponents_complex, exponents_conj)))
    length(exponents_complex) == length(exponents_conj) ||
        throw(ArgumentError("Complex and conjugate exponent lengths are different"))
    length(exponents_complex) == Nc ||
        throw(ArgumentError("Requested $Nc complex variables, but got $(length(exponents_complex))"))
    degree = sum(exponents_complex, init=0) + sum(exponents_conj, init=0)
    idx = exponents_to_index(e, (x[i] for i in 1:Nc for x in (exponents_complex, exponents_conj)), degree)
    @boundscheck iszero(idx) && throw(BoundsError(e, (exponents_complex, exponents_conj)))
    return SimpleMonomial{0,Nc}(unsafe, e, idx, degree)
end

struct OrderedExponents{T,R,C,Cj}
    exponents_real::R
    exponents_complex::C
    exponents_conj::Cj

    OrderedExponents(real::R, complex::C, conj::Cj) where {R<:AbstractVector,C<:AbstractVector,Cj<:AbstractVector} =
        new{promote_type(eltype(R),eltype(C),eltype(Cj)),R,C,Cj}(real, complex, conj)
end

Base.IteratorSize(::Type{<:OrderedExponents}) = Base.HasLength()
Base.IteratorEltype(::Type{<:OrderedExponents}) = Base.HasEltype()
Base.length(oe::OrderedExponents) = length(oe.exponents_real) + 2length(oe.exponents_complex)
Base.eltype(::Type{<:OrderedExponents{T}}) where {T} = T
function Base.iterate(oe::OrderedExponents{T}, (pos, idx)=(0x1, 1)) where {T}
    if pos == 0x1
        if idx > length(oe.exponents_real)
            pos += one(UInt8)
            idx = 1
        else
            @inbounds return oe.exponents_real[idx], (pos, idx +1)
        end
    end
    idx > length(oe.exponents_complex) && return nothing
    if pos == 0x2
        @inbounds return oe.exponents_complex[idx], (0x3, idx)
    else
        @inbounds return oe.exponents_conj[idx], (0x2, idx +1)
    end
end

function SimpleMonomial{Nr,Nc}(e::AbstractExponents{N,I}, exponents_real::AbstractVector{<:Integer},
    exponents_complex::AbstractVector{<:Integer}, exponents_conj::AbstractVector{<:Integer}) where {Nr,Nc,N,I<:Integer}
    N == Nr + 2Nc || throw(MethodError(SimpleMonomial{Nr,Nc}, (e, exponents_real, exponents_complex, exponents_conj)))
    length(exponents_real) == Nr || throw(ArgumentError("Requested $Nr real variables, but got $(length(exponents_real))"))
    length(exponents_complex) == length(exponents_conj) ||
        throw(ArgumentError("Complex and conjugate exponent lengths are different"))
    length(exponents_complex) == Nc ||
        throw(ArgumentError("Requested $Nc complex variables, but got $(length(exponents_complex))"))
    degree = sum(exponents_real, init=0) + sum(exponents_complex, init=0) + sum(exponents_conj, init=0)
    idx = exponents_to_index(e, OrderedExponents(exponents_real, exponents_complex, exponents_conj), degree)
    @boundscheck iszero(idx) && throw(BoundsError(e, (exponents_real, exponents_complex, exponents_conj)))
    return SimpleMonomial{Nr,Nc}(unsafe, e, idx, degree)
end

SimpleMonomial{Nr,Nc}(args::AbstractVector...) where {Nr,Nc} = SimpleMonomial{Nr,Nc,UInt}(args...)
SimpleMonomial{Nr,Nc,I}(args::AbstractVector...) where {Nr,Nc,I<:Integer} =
    SimpleMonomial{Nr,Nc}(ExponentsAll{Nr+2Nc,I}(), args...)

@inline function SimpleMonomial(e::AbstractExponents{N,I}, v::SimpleVariable{Nr,Nc}) where {N,I<:Integer,Nr,Nc}
    N == Nr + 2Nc || throw(MethodError(SimpleMonomial, (e, v)))
    @boundscheck e isa AbstractExponentsDegreeBounded && (e.mindeg > 1 || e.maxdeg < 1) && throw(BoundsError(e, (v,)))
    index_counts(e, 1) # populate the cache
    if e isa ExponentsAll
        return SimpleMonomial{Nr,Nc}(unsafe, e, I(Nr + 2Nc - v.index +2), 1)
    elseif e isa ExponentsDegree
        i = I(Nr + 2Nc - v.index +1)
        iszero(e.mindeg) && (i += one(I))
        return SimpleMonomial{Nr,Nc}(unsafe, e, i, 1)
    else
        @boundscheck @inbounds begin
            (iszero(e.maxmultideg[v.index]) || (e.mindeg > 1) || (isone(e.mindeg) && iszero(e.minmultideg[v.index]))) &&
                throw(BoundsError(e, (v,)))
        end
        isone(e.mindeg) && return SimpleMonomial{Nr,Nc}(unsafe, e, one(I), 1)
        i = I(Nr + 2Nc - v.index +2)
        @inbounds for j in v.index+1:Nr+2Nc
            iszero(e.maxmultideg[j]) && (i -= one(I))
        end
        return SimpleMonomial{Nr,Nc}(unsafe, e, i, 1)
    end
end

"""
    SimpleConjMonomial(m::SimpleMonomial) <: AbstractMonomial

This is a wrapper type for the conjugate of a simple monomial. A lot of operations allow to pass either `SimpleConjMonomial`
or `SimpleMonomial`. Constructing the conjugate using this type works in zero time.

See also [`conj`](@ref Base.conj(::SimpleMonomialOrConj)).
"""
struct SimpleConjMonomial{Nr,Nc,I<:Integer,E<:AbstractExponents} <: AbstractMonomial
    m::SimpleMonomial{Nr,Nc,I,E}

    # don't create conjugates of real monomials
    SimpleConjMonomial(m::SimpleMonomial{Nr,Nc,I,E}) where {Nr,Nc,I<:Integer,E<:AbstractExponents} =
        iszero(Nc) ? m : new{Nr,Nc,I,E}(m)
end

Base.getproperty(c::SimpleConjMonomial, f::Symbol) = getproperty(getfield(c, :m), f)
Base.propertynames(c::SimpleConjMonomial, args...) = propertynames(parent(c), args...)
Base.parent(c::SimpleConjMonomial) = getfield(c, :m)
SimpleConjMonomial(m::SimpleConjMonomial) = parent(m)

"""
    SimpleMonomial(c::SimpleConjMonomial)

Converts a [`SimpleConjMonomial`](@ref) into a [`SimpleMonomial`](@ref). This performs the calculation of the conjugate index.
"""
function SimpleMonomial(c::SimpleConjMonomial{Nr,Nc,<:Integer,<:AbstractExponents}) where {Nr,Nc}
    new_index = exponents_to_index(c.e, exponents(c), degree(c))
    if c.e isa ExponentsMultideg && iszero(new_index) # this can only happen in the multideg case
        throw(ArgumentError("The exponent set does not contain the conjugate monomial"))
    end
    SimpleMonomial{Nr,Nc}(unsafe, c.e, new_index, degree(c))
end
Base.convert(::Type{<:Union{SimpleMonomial,SimpleMonomial{Nr,Nc},SimpleMonomial{Nr,Nc,I},SimpleMonomial{Nr,Nc,I,E}}},
    c::SimpleConjMonomial{Nr,Nc,I,E}) where {Nr,Nc,I<:Integer,E<:AbstractExponents} = SimpleMonomial(c)

MultivariatePolynomials.monomial(m::SimpleMonomial) = m
MultivariatePolynomials.monomial(m::SimpleConjMonomial) = SimpleMonomial(m)

const SimpleMonomialOrConj{Nr,Nc,I<:Integer,E<:AbstractExponents} =
    Union{SimpleMonomial{Nr,Nc,I,E},SimpleConjMonomial{Nr,Nc,I,E}}

Base.isless(x::SimpleMonomial{Nr,Nc}, y::SimpleMonomial{Nr,Nc}) where {Nr,Nc} =
    degree(x) == degree(y) ? compare_indices(unsafe, x.e, x.index, <, y.e, y.index, degree(x)) : isless(degree(x), degree(y))
function Base.isless(x::SimpleMonomialOrConj{Nr,Nc}, y::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc}
    degree(x) == degree(y) || return isless(degree(x), degree(y))
    for (xᵢ, yᵢ) in zip(exponents(x), exponents(y))
        if xᵢ > yᵢ
            return false
        elseif xᵢ < yᵢ
            return true
        end
    end
    return false
end

Base.:(==)(x::SimpleMonomial{Nr,Nc}, y::SimpleMonomial{Nr,Nc}) where {Nr,Nc} =
    degree(x) == degree(y) && compare_indices(unsafe, x.e, x.index, ==, y.e, y.index, degree(x))
function Base.:(==)(x::SimpleMonomialOrConj{Nr,Nc}, y::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc}
    degree(x) == degree(y) || return false
    if (x isa SimpleMonomial && y isa SimpleMonomial) || (x isa SimpleConjMonomial && y isa SimpleConjMonomial)
        return compare_indices(unsafe, x.e, x.index, ==, y.e, y.index, degree(x))
    else
        return all(splat(==), zip(exponents(x), exponents(y)))
    end
end

MultivariatePolynomials.variables(::XorTX{<:SimpleMonomialOrConj{Nr,Nc}}) where {Nr,Nc} = SimpleVariables{Nr,Nc}()
MultivariatePolynomials.nvariables(::XorTX{<:SimpleMonomialOrConj{Nr,Nc}}) where {Nr,Nc} = Nr + 2Nc

function MultivariatePolynomials.convert_constant(T::Type{<:SimpleMonomialOrConj}, α)
    isone(α) || error("Cannot convert `$α` to a `SimpleMonomial` as it is not one")
    return constant_monomial(T)
end

MultivariatePolynomials.degree(m::SimpleMonomialOrConj) = m.degree
MultivariatePolynomials.degree(m::SimpleMonomialOrConj{Nr,Nc}, v::SimpleVariable{Nr,Nc}) where {Nr,Nc} =
    @inbounds exponents(m)[v.index]

# New definition according to https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/pull/292.
for fn in (:degree_complex, :halfdegree)
    @eval function MultivariatePolynomials.$fn(m::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc}
        exps = exponents(parent(m)) # iteration for SimpleMonomial is faster than for conjugate and it doesn't matter here
        Σreal::Int = 0
        Σcomplex::Int = 0
        Σconj::Int = 0
        iter = iterate(exps)::Tuple
        for i in 1:Nr
            Σreal += iter[1]
            if iszero(Nc)
                iter = iterate(exps, iter[2])
                isnothing(iter) && return $(fn === :degree_complex ? :(Σreal) : :(div(Σreal, 2, RoundUp)))
            else
                iter = iterate(exps, iter[2])::Tuple
            end
        end
        for i in 1:Nc
            Σcomplex += iter[1]
            iter = iterate(exps, iter[2])::Tuple
            Σconj += iter[1]
            iter = iterate(exps, iter[2])
            isnothing(iter) && return $(fn === :degree_complex ? :(Σreal) : :(div(Σreal, 2, RoundUp))) + max(Σcomplex, Σconj)
        end
        @assert(false)
    end
end
function MultivariatePolynomials.degree_complex(m::SimpleMonomialOrConj{Nr,Nc}, v::SimpleVariable{Nr,Nc}) where {Nr,Nc}
    exps = exponents(parent(m))
    # while we could use getindex, this will implicitly iterate. So better if we iterate by ourselves, saving some duplication.
    ind₁, ind₂ = minmax(v.index, conj(v).index)
    iter = iterate(exps)::Tuple
    for i in 1:ind₁-1
        iter = iterate(exps, iter[2])::Tuple
    end
    result = iter[1]
    ind₂ == ind₁ && return result
    for i in ind₁:ind₂-1
        iter = iterate(exps, iter[2])::Tuple
    end
    return max(result, iter[1])
end

#region exponents iterator
struct SimpleMonomialExponents{Nr,Nc,Conj,EI<:ExponentIndices} <: AbstractVector{Int}
    ei::EI

    function SimpleMonomialExponents{Nr,Nc}(conj::Bool, ei::EI) where {Nr,Nc,EI<:ExponentIndices}
        length(ei) == Nr + 2Nc || throw(MethodError(SimpleMonomialExponents{Nr,Nc}, (conj, ei)))
        new{Nr,Nc,conj,EI}(ei)
    end
end

Base.IndexStyle(::Type{<:SimpleMonomialExponents}) = IndexLinear()
Base.size(::SimpleMonomialExponents{Nr,Nc}) where {Nr,Nc} = (Nr + 2Nc,)

Base.@propagate_inbounds function Base.getindex(sme::SimpleMonomialExponents{Nr,Nc,Conj}, varidx::Integer) where {Nr,Nc,Conj}
    if Conj
        if varidx > Nr
            varidx = Nr + one(Nr) + ((varidx - Nr - one(Nr)) ⊻ one(varidx))
        end
    end
    return sme.ei[varidx]
end

Base.iterate(sme::SimpleMonomialExponents{<:Any,<:Any,false}, args...) = iterate(sme.ei, args...)
function Base.iterate(sme::SimpleMonomialExponents{Nr,Nc,true}) where {Nr,Nc}
    @assert(!iszero(Nc))
    next = iterate(sme.ei)::Tuple
    # We want to keep the internal interface that the first entry in the iterator state is the remaining degree, so a little
    # bit of unpacking is necessary.
    if iszero(Nr)
        second = iterate(sme.ei, next[2])::Tuple
        return second[1], (second[2]..., next[1])
    else
        return next[1], (next[2]..., -1)
    end
end
function Base.iterate(sme::SimpleMonomialExponents{Nr,Nc,true}, allstate) where {Nr,Nc}
    state = allstate[1:end-1]
    prev = last(allstate)
    # we know that Nc ≠ 0, for it is not possible to create a real-valued SimpleConjMonomial
    @assert(!iszero(Nc))
    i = state[2] # standardized part of state
    if i ≤ Nr +1
        next = iterate(sme.ei, state)::Tuple
        return next[1], (next[2]..., -1)
    elseif prev == -1
        next = iterate(sme.ei, state)
        isnothing(next) && return nothing # might happen when we are at the end
        second = iterate(sme.ei, next[2])::Tuple # but then this cannot happen
        return second[1], (second[2]..., next[1])
    else
        # everything was already calculated previously, just yield it and go on
        return prev, (state..., -1)
    end
end
# iteration is faster than indexed access, so let's fall back to the generic iterator-based functions
Base.@propagate_inbounds Base.copyto!(dest::AbstractArray, src::SimpleMonomialExponents) =
    @invoke copyto!(dest::AbstractArray, src::Any)
Base.@propagate_inbounds Base.copyto!(dest::AbstractArray, dstart::Integer, src::SimpleMonomialExponents) =
    @invoke copyto!(dest::AbstractArray, dstart::Integer, src::Any)
Base.@propagate_inbounds Base.copyto!(dest::AbstractArray, dstart::Integer, src::SimpleMonomialExponents, sstart::Integer) =
    @invoke copyto!(dest::AbstractArray, dstart::Integer, src::Any, sstart::Integer)
Base.@propagate_inbounds Base.copyto!(dest::AbstractArray, dstart::Integer, src::SimpleMonomialExponents, sstart::Integer,
    n::Integer) = @invoke copyto!(dest::AbstractArray, dstart::Integer, src::Any, sstart::Integer, n::Integer)

Base.sum(sme::SimpleMonomialExponents; init=0) = sum(sme.ei; init)

MultivariatePolynomials._zip(t::Tuple, e::SimpleMonomialExponents) = zip(t, e)
#endregion
MultivariatePolynomials.exponents(m::SimpleMonomialOrConj{Nr,Nc,I}) where {Nr,Nc,I<:Integer} =
    SimpleMonomialExponents{Nr,Nc}(
        m isa SimpleConjMonomial,
        exponents_from_index(unsafe, m.e, m.index, degree(m))
    )

# implement an iteration method although there is the exponents function - this one gives a (SimpleVariable, exponent) tuple
# and skips over zero exponents
Base.IteratorSize(::Type{<:SimpleMonomialOrConj}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:SimpleMonomialOrConj}) = Base.HasEltype()
Base.eltype(::Type{<:SimpleMonomialOrConj{Nr,Nc}}) where {Nr,Nc} =
    Tuple{SimpleVariable{Nr,Nc,smallest_unsigned(Nr + 2Nc)},Int}
function Base.iterate(m::SimpleMonomialOrConj{Nr,Nc}) where {Nr,Nc}
    sme = exponents(m)
    iter = iterate(sme)::Tuple
    exponent, state = iter
    if iszero(exponent)
        return iterate(m, (sme, 1, state))
    else
        return (SimpleVariable{Nr,Nc}(1), exponent), (sme, 1, state)
    end
end
function Base.iterate(m::SimpleMonomialOrConj{Nr,Nc}, (sme, lastvar, state)) where {Nr,Nc}
    # state[1]: remaining degree
    while true
        iszero(state[1]) && (!(m isa SimpleConjMonomial) || last(state) < 1) && return nothing
        # ^ fast-path: if the remaining degree is already zero, we can finish; but note that for the conj monomial, this will
        # be the remaining degree after the current ordinary variable and its conjugate have been traveled.
        iter = iterate(sme, state)
        isnothing(iter) && return nothing
        lastvar += 1
        exponent, state = iter
        iszero(exponent) || return (SimpleVariable{Nr,Nc}(lastvar), exponent), (sme, lastvar, state)
        # index of next variable after iteration -> current variable is -2
    end
end

MultivariatePolynomials.isconstant(m::SimpleMonomialOrConj{<:Any,<:Any,<:Integer,<:ExponentsAll}) = isone(m.index)
MultivariatePolynomials.isconstant(m::SimpleMonomialOrConj{<:Any,<:Any,<:Integer,<:AbstractExponentsDegreeBounded}) =
    isone(m.index) && iszero(m.e.mindeg)

MultivariatePolynomials.constant_monomial(m::SimpleMonomialOrConj{Nr,Nc,I,<:ExponentsAll}) where {Nr,Nc,I<:Integer} =
    SimpleMonomial{Nr,Nc}(unsafe, m.e, one(I), 0)
MultivariatePolynomials.constant_monomial(m::SimpleMonomialOrConj{Nr,Nc,I,<:AbstractExponentsDegreeBounded}) where {Nr,Nc,I<:Integer} =
    iszero(m.e.mindeg) ? SimpleMonomial{Nr,Nc}(unsafe, e, one(I), 0) :
                         throw(ArgumentError("Constant monomial is not part of the exponent set"))
function MultivariatePolynomials.constant_monomial(::Type{<:Union{<:SimpleMonomialOrConj{Nr,Nc,I_},
                                                                  <:SimpleMonomialOrConj{Nr,Nc}}}) where {Nr,Nc,I_<:Integer}
    I = @isdefined(I_) ? I_ : UInt
    SimpleMonomial{Nr,Nc}(unsafe, ExponentsAll{Nr+2Nc,I}(), one(I), 0) # we cannot obtain the necessary information about the
                                                                       # exponents, so we must change the type
end

Base.conj(m::SimpleMonomial{Nr,0} where {Nr}) = m
Base.conj(m::SimpleMonomial{Nr,Nc}) where {Nr,Nc} = SimpleMonomial(SimpleConjMonomial(m))
Base.conj(m::SimpleConjMonomial{Nr,Nc}) where {Nr,Nc} = parent(m)
"""
    conj(m::Union{<:SimpleMonomial,<:SimpleConjMonomial})

Creates the conjugate of a [`SimpleMonomial`](@ref). The result type of this operation will always be [`SimpleMonomial`](@ref).
If the conjugate can be used to work with lazily, consider wrapping the monomial in a [`SimpleConjMonomial`](@ref) instead.
"""
Base.conj(::SimpleMonomialOrConj)

Base.@assume_effects :consistent function Base.isreal(m::SimpleMonomial{Nr,Nc}) where {Nr,Nc}
    exps = exponents(m)
    iter = iterate(exps)
    for _ in 1:Nr
        isnothing(iter) && return true # should not happen, but let's help the compiler to figure out type stability
        iszero(iter[2][1]) && return true # short-circuit if no nonzero exponents follow
        iter = iterate(exps, iter[2])
    end
    for _ in 1:Nc
        isnothing(iter) && return true
        complexdeg = iter[1]
        iter = iterate(exps, iter[2])::Tuple
        iter[1] == complexdeg || return false
        iszero(iter[2][1]) && return true # short-circuit if the remaining degree is zero
        isodd(iter[2][1]) && return false # short-circuit if an odd number of complex variables remain, they can never be
                                          # distributed evenly.
        iter = iterate(exps, iter[2])
    end
    return true
end
Base.isreal(m::SimpleConjMonomial) = isreal(parent(m)) # iteration for the ordinary monomial is faster

# we hash the monomial to its index, assuming that usually we are interested in a list of monomials of the same type, so no
# collisions should occur in this way.
Base.hash(m::SimpleMonomialOrConj, h::UInt) = hash(m.index, h)

effective_variables_in(m::SimpleMonomial, in) = all(vp -> ordinary_variable(vp[1]) ∈ in, m)

"""
    monomial_product(e::AbstractExponents, m...)

Calculates the product of all monomials (or conjugates, or variables) `m`. The result must be part of the exponent set `e`.
If the default multiplication `*` is used instead, `e` will always be `ExponentsAll` with the jointly promoted index type.
"""
Base.@assume_effects :consistent function monomial_product(e::AbstractExponents{N},
                                                           m::Union{<:SimpleMonomialOrConj{Nr,Nc},
                                                                    <:SimpleVariable{Nr,Nc}}...) where {N,Nr,Nc}
    N == Nr + 2Nc || throw(MethodError(monomial_product, (e, m...)))
    index, d = exponents_sum(e, exponents.(m)...)
    e isa ExponentsAll ||
        (iszero(index) && throw(ArgumentError("The given product is not present in the required exponent range")))
    return SimpleMonomial{Nr,Nc}(unsafe, e, index, d)
end

_get_I(::Type{<:SimpleMonomialOrConj{<:Any,<:Any,I}}) where {I<:Integer} = I
_get_I(::Type{<:SimpleVariable}) = missing
@generated function _get_I(m::Union{<:SimpleMonomialOrConj{Nr,Nc},<:SimpleVariable{Nr,Nc}}...) where {Nr,Nc}
    I = missing
    for T in m
        if ismissing(I)
            I = _get_I(T)
        else
            newI = _get_I(T)
            if !ismissing(newI)
                I = promote_type(I, newI)
            end
        end
    end
    if ismissing(I)
        I = UInt
    end
    return :(return $I)
end

Base.:*(m::Union{<:SimpleMonomialOrConj{Nr,Nc},<:SimpleVariable{Nr,Nc}}...) where {Nr,Nc} =
    monomial_product(ExponentsAll{Nr+2Nc,_get_I(m...)}(), m...)

function Base.:^(m::SimpleMonomialOrConj{Nr,Nc}, p::Integer) where {Nr,Nc}
    idx, deg = exponents_product(m.e, exponents(m), p)
    @boundscheck !(m.e isa ExponentsAll) && iszero(idx) && throw(BoundsError(m.e, idx))
    SimpleMonomial{Nr,Nc}(unsafe, m.e, idx, deg)
end

"""
    monomial_index([e::AbstractExponents,] m...)

Calculates the index of the given monomial (or the product of all given monomials, or conjugates, or variables) `m`. The result
must be part of the exponent set `e`. If `e` is omitted, it will be be `ExponentsAll` with the jointly promoted index type.
"""
@inline function monomial_index(e::AbstractExponents{N}, m::SimpleMonomial{<:Any,<:Any,<:Integer,<:AbstractExponents{N}}) where {N}
    e === m.e && return m.index
    d = degree(m)
    index_counts(e, d) # we don't know whether the cache in `e` is sufficient
    return convert_index(unsafe, e, m.e, m.index, d)
end

@inline function monomial_index(e::AbstractExponents{N}, m::SimpleConjMonomial{<:Any,<:Any,<:Integer,<:AbstractExponents{N}}) where {N}
    d = degree(m)
    e === m.e || index_counts(e, d)
    return exponents_to_index(e, exponents(m), d)
end

monomial_index(m::Union{<:SimpleMonomialOrConj{Nr,Nc},<:SimpleVariable{Nr,Nc}}...) where {Nr,Nc} =
    monomial_index(ExponentsAll{Nr+2Nc,_get_I(m...)}(), m...)
monomial_index(e::AbstractExponents, m::Union{<:SimpleMonomialOrConj{Nr,Nc},<:SimpleVariable{Nr,Nc}}...) where {Nr,Nc} =
    monomial_product(e, m...).index
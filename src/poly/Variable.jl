export IntVariable, IntRealVariable, IntComplexVariable, variable_index

"""
    IntVariable{Nr,Nc,I<:Unsigned} <: AbstractVariable

`IntVariable` is the basic type for any simple variable in a polynomial right with `Nr` real and `Nc` complex-valued
variables. The variable is identified by its index of type `I` alone. A variable can be explicitly cast to `I` in order to
obtain its index.

To construct a real-valued or complex-valued variable, see [`IntRealVariable`](@ref) and [`IntComplexVariable`](@ref).

!!! warning
    Note that `I` has nothing to do with the index type used to identify monomials or exponents; in fact, `I` is automatically
    calculated as the smallest `Unsigned` descendant that can still hold the value `Nr+2Nc`.
"""
struct IntVariable{Nr,Nc,I<:Unsigned} <: AbstractVariable
    index::I

    @inline function IntVariable{Nr,Nc}(index::Integer) where {Nr,Nc}
        0 ≤ Nr || throw(DomainError(Nr, "Invalid ring: Nr must be a nonnegative integer"))
        0 ≤ Nc || throw(DomainError(Nc, "Invalid ring: Nc must be nonnegative integer"))
        @boundscheck 0 < index ≤ Nr + 2Nc || throw(DomainError(index, "Invalid index: must be between 1 and $(Nr + 2Nc)"))
        I = smallest_unsigned(Nr + 2Nc)
        new{Nr,Nc,I}(I(index))
    end
end

(::Type{I})(v::IntVariable{<:Any,<:Any,I}) where {I<:Unsigned} = v.index

# we want curly-initialization style, so use this singleton that can never be constructed
struct IntRealVariable{Nr,Nc}
    @doc """
        IntRealVariable{Nr,Nc}(index::Integer)

    Creates a new real-valued simple variable whose identity is determined by `index`. Real-valued variables with the same
    index are considered identical; however, they are different from complex-valued variables constructed with the same index.
    A real variable will print as `xᵢ`, where the subscript is given by `index`.
    The variable is part of the polynomial ring with `Nr` real and `Nc` complex variables.

    !!! warning
        This method is for construction of the variable only. Do not use it in type comparisons; variables constructed with
        this method will _not_ be of type `IntRealVariable` (in fact, don't think of it as a type), but rather of type
        [`IntVariable`](@ref)!

    See also [`IntVariable`](@ref), [`IntComplexVariable`](@ref).
    """
    function IntRealVariable{Nr,Nc}(index::Integer) where {Nr,Nc}
        0 < index ≤ Nr || throw(DomainError(index, "Invalid index: must be between 1 and $Nr"))
        return IntVariable{Nr,Nc}(index)
    end
end

struct IntComplexVariable{Nr,Nc}
    @doc """
        IntComplexVariable{Nr,Nc}(index::Integer, isconj::Bool=false)

    Creates a new complex-valued simple variable whose identity is determined by `index`, and which is a conjugate variable if
    `isconj` is set appropriately. Complex-valued variables with the same index and `isconj` state are considered identical;
    however, they are different from real-valued variables constructed with the same index, even if they are not conjugate.
    A complex variable will print as `zᵢ` (if `isconj=false`) or `z̄ᵢ` (if `isconj=true`), where the subscript is given by
    `index`.
    The variable is part of the polynomial ring with `Nr` real and `Nc` complex variables.

    !!! warning
        This method is for construction of the variable only. Do not use it in type comparisons; variables constructed with
        this method will _not_ be of type `IntComplexVariable` (in fact, don't think of it as a type), but rather of type
        [`IntVariable`](@ref)!

    See also [`IntVariable`](@ref), [`IntRealVariable`](@ref).
    """
    function IntComplexVariable{Nr,Nc}(index::Integer, isconj::Bool=false) where {Nr,Nc}
        0 < index ≤ Nc || throw(DomainError(index, "Invalid index: must be between 1 and $Nc"))
        return IntVariable{Nr,Nc}(2index + Nr - typeof(index)(!isconj))
    end
end

function MultivariatePolynomials.name(v::IntVariable{Nr,Nc}) where {Nr,Nc}
    if iszero(Nc) || v.index ≤ Nr # compile-time short-circuit
        return "x" * MultivariatePolynomials.unicode_subscript(v.index)
    else
        cindex = v.index - Nr
        return "z" * MultivariatePolynomials.unicode_subscript(div(cindex, 2, RoundUp), iseven(cindex))
    end
end

function MultivariatePolynomials.name_base_indices(v::IntVariable{Nr,Nc}) where {Nr,Nc}
    if iszero(Nc) || v.index ≤ Nr # compile-time short-circuit
        return (:x, v.index)
    else
        cindex = v.index - Nr
        return (:z, div(cindex, 2, RoundUp), iseven(cindex))
    end
end

Base.:(==)(x::V, y::V) where {V<:IntVariable} = x.index == y.index
Base.isless(x::V, y::V) where {V<:IntVariable} = isless(x.index, y.index)

Base.IteratorSize(::Type{<:IntVariable}) = Base.HasLength()
Base.IteratorEltype(::Type{<:IntVariable}) = Base.HasEltype()
Base.eltype(::Type{V}) where {V<:IntVariable} = Tuple{V,Int}
Base.length(::IntVariable) = 1
Base.iterate(v::IntVariable) = (v, 1), missing
Base.iterate(::IntVariable, ::Missing) = nothing

#region Variables iterator
struct IntVariables{Nr,Nc,V<:IntVariable{Nr,Nc}} <: AbstractVector{V}
    IntVariables{Nr,Nc}() where {Nr,Nc} = new{Nr,Nc,IntVariable{Nr,Nc,smallest_unsigned(Nr + 2Nc)}}()
end

Base.IndexStyle(::Type{<:IntVariables}) = IndexLinear()
Base.size(::IntVariables{Nr,Nc}) where {Nr,Nc} = (Nr + 2Nc,)

@inline function Base.getindex(smv::IntVariables{Nr,Nc}, idx::Integer) where {Nr,Nc}
    @boundscheck checkbounds(smv, idx)
    return IntVariable{Nr,Nc}(idx)
end

Base.collect(::IntVariables{Nr,Nc}) where {Nr,Nc} = map(IntVariable{Nr,Nc}, 1:Nr+2Nc)
#endregion

MultivariatePolynomials.variables(::XorTX{<:IntVariable{Nr,Nc}}) where {Nr,Nc} = IntVariables{Nr,Nc}()
MultivariatePolynomials.nvariables(::XorTX{<:IntVariable{Nr,Nc}}) where {Nr,Nc} = Nr + 2Nc

function MultivariatePolynomials.monomial(v::IntVariable{Nr,Nc}) where {Nr,Nc}
    # construct the monomial directly. But beware that the constructor is unsafe, so we need to make sure the cache is set up
    # properly before. Also, we don't know the exponents type, so this function just uses a default.
    e = ExponentsAll{Nr+2Nc,UInt}()
    index_counts(e, 1) # populate cache
    return IntMonomial{Nr,Nc}(unsafe, e, UInt(Nr + 2Nc - v.index +2), 1)
end

#region exponents iterator
struct IntVariableExponents{Nr,Nc,V<:IntVariable{Nr,Nc}} <: AbstractVector{Int}
    v::V
end

Base.IndexStyle(::Type{<:IntVariableExponents}) = IndexLinear()
Base.size(::IntVariableExponents{Nr,Nc}) where {Nr,Nc} = (Nr + 2Nc,)

@inline function Base.getindex(sve::IntVariableExponents, idx::Integer)
    @boundscheck checkbounds(sve, idx)
    return Int(idx == sve.v.index)
end

function Base.collect(sve::IntVariableExponents{Nr,Nc}) where {Nr,Nc}
    result = zeros(Int, Nr + 2Nc)
    @inbounds result[sve.v.index] = 1
    return result
end

MultivariatePolynomials._zip(t::Tuple, e::IntVariableExponents) = zip(t, e)
#endregion
MultivariatePolynomials.exponents(v::IntVariable) = IntVariableExponents(v)

Base.conj(v::IntVariable{<:Any,0}) = v
Base.conj(v::IntVariable{Nr,Nc}) where {Nr,Nc} =
    IntVariable{Nr,Nc}(v.index ≤ Nr ? v.index : (Nr + one(Nr) + ((v.index - Nr - one(Nr)) ⊻ one(v.index))))

Base.real(::IntVariable) = error("Not implemented")
Base.imag(::IntVariable) = error("Not implemented")

MultivariatePolynomials.isreal(::IntVariable{<:Any,0}) = true
MultivariatePolynomials.isreal(v::IntVariable{Nr}) where {Nr} = v.index ≤ Nr
MultivariatePolynomials.isconj(::IntVariable{<:Any,0}) = false
MultivariatePolynomials.isconj(v::IntVariable{Nr,Nc}) where {Nr,Nc} = v.index > Nr && iseven(v.index - Nr)

MultivariatePolynomials.ordinary_variable(v::IntVariable{<:Any,0}) = v
MultivariatePolynomials.ordinary_variable(v::IntVariable{Nr,Nc}) where {Nr,Nc} =
    IntVariable{Nr,Nc}(v.index ≤ Nr ? v.index : (Nr + ((v.index - Nr - one(Nr)) | one(v.index))))

"""
    variable_index(v::IntVariable{Nr,Nc})

Returns the index of the variable `v`, where real-valued variables have indices between 1 and `Nr`, and complex-valued
variables have indices between `Nr+1` and `Nr+Nc`. Conjugates have the same indices as their ordinary variables.
"""
variable_index(v::IntVariable{<:Any,0}) = v.index
variable_index(v::IntVariable{Nr,<:Any,I}) where {Nr,I<:Integer} =
    v.index ≤ Nr ? v.index : I(Nr + ((v.index - Nr -1) >> 1) + one(I))

function monomial_index(::ExponentsAll{N,I}, v::IntVariable{Nr,Nc}) where {N,I<:Integer,Nr,Nc}
    N == Nr + 2Nc || throw(MethodError(monomial_index, (e, v)))
    return I(N - v.index +2)
end
monomial_index(e::AbstractExponents, v::IntVariable) = exponents_to_index(e, exponents(v), 1)

function Base.:^(v::IntVariable{Nr,Nc}, p::Integer) where {Nr,Nc}
    e = ExponentsAll{Nr+2Nc,UInt}()
    return IntMonomial{Nr,Nc}(unsafe, e, exponents_product(e, exponents(v), p)...)
end
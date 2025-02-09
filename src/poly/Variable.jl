export SimpleVariable, SimpleRealVariable, SimpleComplexVariable, variable_index

"""
    SimpleVariable{Nr,Nc,I<:Unsigned} <: AbstractVariable

`SimpleVariable` is the basic type for any simple variable in a polynomial right with `Nr` real and `Nc` complex-valued
variables. The variable is identified by its index of type `I` alone. A variable can be explicitly cast to `I` in order to
obtain its index.

To construct a real-valued or complex-valued variable, see [`SimpleRealVariable`](@ref) and [`SimpleComplexVariable`](@ref).

!!! warning
    Note that `I` has nothing to do with the index type used to identify monomials or exponents; in fact, `I` is automatically
    calculated as the smallest `Unsigned` descendant that can still hold the value `Nr+2Nc`.
"""
struct SimpleVariable{Nr,Nc,I<:Unsigned} <: AbstractVariable
    index::I

    @inline function SimpleVariable{Nr,Nc}(index::Integer) where {Nr,Nc}
        0 ≤ Nr || throw(DomainError(Nr, "Invalid ring: Nr must be a nonnegative integer"))
        0 ≤ Nc || throw(DomainError(Nc, "Invalid ring: Nc must be nonnegative integer"))
        @boundscheck 0 < index ≤ Nr + 2Nc || throw(DomainError(index, "Invalid index: must be between 1 and $(Nr + 2Nc)"))
        I = smallest_unsigned(Nr + 2Nc)
        new{Nr,Nc,I}(I(index))
    end
end

(::Type{I})(v::SimpleVariable{<:Any,<:Any,I}) where {I<:Unsigned} = v.index

# we want curly-initialization style, so use this singleton that can never be constructed
struct SimpleRealVariable{Nr,Nc}
    @doc """
        SimpleRealVariable{Nr,Nc}(index::Integer)

    Creates a new real-valued simple variable whose identity is determined by `index`. Real-valued variables with the same
    index are considered identical; however, they are different from complex-valued variables constructed with the same index.
    A real variable will print as `xᵢ`, where the subscript is given by `index`.
    The variable is part of the polynomial ring with `Nr` real and `Nc` complex variables.

    !!! warning
        This method is for construction of the variable only. Do not use it in type comparisons; variables constructed with
        this method will _not_ be of type `SimpleRealVariable` (in fact, don't think of it as a type), but rather of type
        [`SimpleVariable`](@ref)!

    See also [`SimpleVariable`](@ref), [`SimpleComplexVariable`](@ref).
    """
    function SimpleRealVariable{Nr,Nc}(index::Integer) where {Nr,Nc}
        0 < index ≤ Nr || throw(DomainError(index, "Invalid index: must be between 1 and $Nr"))
        return SimpleVariable{Nr,Nc}(index)
    end
end

struct SimpleComplexVariable{Nr,Nc}
    @doc """
        SimpleComplexVariable{Nr,Nc}(index::Integer, isconj::Bool=false)

    Creates a new complex-valued simple variable whose identity is determined by `index`, and which is a conjugate variable if
    `isconj` is set appropriately. Complex-valued variables with the same index and `isconj` state are considered identical;
    however, they are different from real-valued variables constructed with the same index, even if they are not conjugate.
    A complex variable will print as `zᵢ` (if `isconj=false`) or `z̄ᵢ` (if `isconj=true`), where the subscript is given by
    `index`.
    The variable is part of the polynomial ring with `Nr` real and `Nc` complex variables.

    !!! warning
        This method is for construction of the variable only. Do not use it in type comparisons; variables constructed with
        this method will _not_ be of type `SimpleComplexVariable` (in fact, don't think of it as a type), but rather of type
        [`SimpleVariable`](@ref)!

    See also [`SimpleVariable`](@ref), [`SimpleRealVariable`](@ref).
    """
    function SimpleComplexVariable{Nr,Nc}(index::Integer, isconj::Bool=false) where {Nr,Nc}
        0 < index ≤ Nc || throw(DomainError(index, "Invalid index: must be between 1 and $Nc"))
        return SimpleVariable{Nr,Nc}(2index + Nr - typeof(index)(!isconj))
    end
end

function MultivariatePolynomials.name(v::SimpleVariable{Nr,Nc}) where {Nr,Nc}
    if iszero(Nc) || v.index ≤ Nr # compile-time short-circuit
        return "x" * MultivariatePolynomials.unicode_subscript(v.index)
    else
        cindex = v.index - Nr
        return "z" * MultivariatePolynomials.unicode_subscript(div(cindex, 2, RoundUp), iseven(cindex))
    end
end

function MultivariatePolynomials.name_base_indices(v::SimpleVariable{Nr,Nc}) where {Nr,Nc}
    if iszero(Nc) || v.index ≤ Nr # compile-time short-circuit
        return (:x, v.index)
    else
        cindex = v.index - Nr
        return (:z, div(cindex, 2, RoundUp), iseven(cindex))
    end
end

Base.:(==)(x::V, y::V) where {V<:SimpleVariable} = x.index == y.index
Base.isless(x::V, y::V) where {V<:SimpleVariable} = isless(x.index, y.index)

#region Variables iterator
struct SimpleVariables{Nr,Nc,V<:SimpleVariable{Nr,Nc}} <: AbstractVector{V}
    SimpleVariables{Nr,Nc}() where {Nr,Nc} = new{Nr,Nc,SimpleVariable{Nr,Nc,smallest_unsigned(Nr + 2Nc)}}()
end

Base.IndexStyle(::Type{<:SimpleVariables}) = IndexLinear()
Base.size(::SimpleVariables{Nr,Nc}) where {Nr,Nc} = (Nr + 2Nc,)

@inline function Base.getindex(smv::SimpleVariables{Nr,Nc}, idx::Integer) where {Nr,Nc}
    @boundscheck checkbounds(smv, idx)
    return SimpleVariable{Nr,Nc}(idx)
end

Base.collect(::SimpleVariables{Nr,Nc}) where {Nr,Nc} = map(SimpleVariable{Nr,Nc}, 1:Nr+2Nc)
#endregion

MultivariatePolynomials.variables(::XorTX{<:SimpleVariable{Nr,Nc}}) where {Nr,Nc} = SimpleVariables{Nr,Nc}()
MultivariatePolynomials.nvariables(::XorTX{<:SimpleVariable{Nr,Nc}}) where {Nr,Nc} = Nr + 2Nc

function MultivariatePolynomials.monomial(v::SimpleVariable{Nr,Nc}) where {Nr,Nc}
    # construct the monomial directly. But beware that the constructor is unsafe, so we need to make sure the cache is set up
    # properly before. Also, we don't know the exponents type, so this function just uses a default.
    e = ExponentsAll{Nr+2Nc,UInt}()
    index_counts(e, 1) # populate cache
    return SimpleMonomial{Nr,Nc}(unsafe, e, UInt(Nr + 2Nc - v.index +2), 1)
end

#region exponents iterator
struct SimpleVariableExponents{Nr,Nc,V<:SimpleVariable{Nr,Nc}} <: AbstractVector{Int}
    v::V
end

Base.IndexStyle(::Type{<:SimpleVariableExponents}) = IndexLinear()
Base.size(::SimpleVariableExponents{Nr,Nc}) where {Nr,Nc} = (Nr + 2Nc,)

@inline function Base.getindex(sve::SimpleVariableExponents, idx::Integer)
    @boundscheck checkbounds(sve, idx)
    return Int(idx == sve.v.index)
end

function Base.collect(sve::SimpleVariableExponents{Nr,Nc}) where {Nr,Nc}
    result = zeros(Int, Nr + 2Nc)
    @inbounds result[sve.v.index] = 1
    return result
end

MultivariatePolynomials._zip(t::Tuple, e::SimpleVariableExponents) = zip(t, e)
#endregion
MultivariatePolynomials.exponents(v::SimpleVariable) = SimpleVariableExponents(v)

Base.conj(v::SimpleVariable{<:Any,0}) = v
Base.conj(v::SimpleVariable{Nr,Nc}) where {Nr,Nc} =
    SimpleVariable{Nr,Nc}(v.index ≤ Nr ? v.index : (Nr + one(Nr) + ((v.index - Nr - one(Nr)) ⊻ one(v.index))))

Base.real(::SimpleVariable) = error("Not implemented")
Base.imag(::SimpleVariable) = error("Not implemented")

MultivariatePolynomials.isreal(::SimpleVariable{<:Any,0}) = true
MultivariatePolynomials.isreal(v::SimpleVariable{Nr}) where {Nr} = v.index ≤ Nr
MultivariatePolynomials.isconj(::SimpleVariable{<:Any,0}) = false
MultivariatePolynomials.isconj(v::SimpleVariable{Nr,Nc}) where {Nr,Nc} = v.index > Nr && iseven(v.index - Nr)

MultivariatePolynomials.ordinary_variable(v::SimpleVariable{<:Any,0}) = v
MultivariatePolynomials.ordinary_variable(v::SimpleVariable{Nr,Nc}) where {Nr,Nc} =
    SimpleVariable{Nr,Nc}(v.index ≤ Nr ? v.index : (Nr + ((v.index - Nr - one(Nr)) | one(v.index))))

"""
    variable_index(v::SimpleVariable{Nr,Nc})

Returns the index of the variable `v`, where real-valued variables have indices between 1 and `Nr`, and complex-valued
variables have indices between `Nr+1` and `Nr+Nc`. Conjugates have the same indices as their ordinary variables.
"""
variable_index(v::SimpleVariable{<:Any,0}) = v.index
variable_index(v::SimpleVariable{Nr,<:Any,I}) where {Nr,I<:Integer} =
    v.index ≤ Nr ? v.index : I(Nr + ((v.index - Nr -1) >> 1) + one(I))

function monomial_index(::ExponentsAll{N,I}, v::SimpleVariable{Nr,Nc}) where {N,I<:Integer,Nr,Nc}
    N == Nr + 2Nc || throw(MethodError(monomial_index, (e, v)))
    return I(N - v.index +2)
end
monomial_index(e::AbstractExponents, v::SimpleVariable) = exponents_to_index(e, exponents(v), 1)

function Base.:^(v::SimpleVariable{Nr,Nc}, p::Integer) where {Nr,Nc}
    e = ExponentsAll{Nr+2Nc,UInt}()
    return SimpleMonomial{Nr,Nc}(unsafe, e, exponents_product(e, exponents(v), p)...)
end
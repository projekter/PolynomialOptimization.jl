export SimplePolynomial

struct SimplePolynomial{C,Nr,Nc,M<:SimpleMonomialVector{Nr,Nc}} <: AbstractPolynomial{C}
    coeffs::Vector{C}
    monomials::M

    """
        SimplePolynomial(coeffs::Vector{C}, monomials::SimpleMonomialVector)

    Creates a `SimplePolynomial` with well-defined coefficients and monomials, which must have the same length.
    """
    @inline function SimplePolynomial(coeffs::Vector{C}, monomials::M) where {C,Nr,Nc,M<:SimpleMonomialVector{Nr,Nc}}
        length(coeffs) == length(monomials) || throw(ArgumentError("Lengths are different"))
        return new{C,Nr,Nc,M}(coeffs, monomials)
    end
end

MultivariatePolynomials.variables(p::XorTX{SimplePolynomial}) = variables(monomial_type(p))

MultivariatePolynomials.coefficients(p::SimplePolynomial) = p.coeffs

function MultivariatePolynomials.coefficient(p::SimplePolynomial{C,Nr,Nc}, m::SimpleMonomial{Nr,Nc}) where {C,Nr,Nc}
    pos = searchsortedlast(p.monomials, m)
    (iszero(pos) || @inbounds(p.monomials[pos] != m)) && return zero(C)
    @inbounds return p.coeffs[pos]
end

MultivariatePolynomials.monomials(p::SimplePolynomial) = p.monomials

MultivariatePolynomials.terms(p::SimplePolynomial) = collect(p)

MultivariatePolynomials.nterms(p::SimplePolynomial) = length(p.coeffs)

MultivariatePolynomials.nvariables(::XorTX{SimplePolynomial{C,Nr,Nc} where {C}}) where {Nr,Nc} = Nr + 2Nc

Base.iterate(p::SimplePolynomial) =
    isempty(p.coeffs) ? nothing : (@inbounds(term_type(p)(p.coeffs[begin], p.monomials[begin])), 2)
Base.iterate(p::SimplePolynomial, i::Integer) =
    i > length(p.coeffs) ? nothing : (@inbounds(term_type(p)(p.coeffs[i], p.monomials[i])), i +1)
Base.IteratorSize(::Type{<:SimplePolynomial}) = Base.HasLength()
Base.IteratorEltype(::Type{<:SimplePolynomial}) = Base.HasEltype()
Base.eltype(::Type{P}) where {P<:SimplePolynomial} = term_type(P)
Base.length(p::SimplePolynomial) = length(p.coeffs)

for f in (:mindegree, :maxdegree, :extdegree, :mindegree_complex, :maxdegree_complex, :extdegree_complex, :minhalfdegree,
          :maxhalfdegree, :exthalfdegree, :variables)
    @eval MultivariatePolynomials.$f(p::SimplePolynomial) = $f(p.monomials)
end

MultivariatePolynomials.constant_monomial(m::SimplePolynomial) = constant_monomial(first(m.monomials))

MultivariatePolynomials.constant_monomial(P::Type{<:SimplePolynomial}) = constant_monomial(eltype(P))

MultivariatePolynomials.map_coefficients(f::Function, p::SimplePolynomial) = SimplePolynomial(map(f, p.coeffs), p.monomials)

function MultivariatePolynomials.map_coefficients!(f::Function, p::SimplePolynomial)
    @inbounds for i in 1:length(p.coeffs)
        p.coeffs[i] = f(p.coeffs[i])
    end
    return p
end

function MultivariatePolynomials.effective_variables(p::SimplePolynomial{<:Any,Nr,Nc}; rettype::Type{V}=Vector, by::F=identity) where {Nr,Nc,V<:Union{Vector,Set},F}
    ev = Set{Base.promote_op(by, SimpleVariable{Nr,Nc,smallest_unsigned(Nr+2Nc)})}()
    for t in p, (v, _) in monomial(t)
        push!(ev, by(v))
    end
    V <: Vector && return sort!(collect(ev))
    return ev
end

# these are tricky. They shouldn't even be there since they might (or will) violate the no-allocation policy of
# SimplePolynomials
function Base.conj(p::SimplePolynomial)
    newcoeffs = conj(p.coeffs)
    SimplePolynomial(newcoeffs, conj(p.monomials, newcoeffs)) # conjugation implies reordering, so sort the coefficients along
end
function MultivariatePolynomials.LinearAlgebra.adjoint(p::SimplePolynomial)
    newcoeffs = adjoint.(p.coeffs)
    SimplePolynomial(newcoeffs, conj(p.monomials, newcoeffs))
end
Base.isreal(p::SimplePolynomial{<:Any,<:Any,0}) = all(∘(iszero, imag), p.coeffs)
function Base.isreal(p::SimplePolynomial)
    # the monomial vector is sorted and unique
    for t in p
        coefficient(t) == conj(coefficient(p, SimpleConjMonomial(monomial(t)))) || return false
    end
    return true
end

MultivariatePolynomials.LinearAlgebra.transpose(p::SimplePolynomial) = SimplePolynomial(transpose.(p.coeffs), p.monomials)

SimplePolynomial(p::SimplePolynomial) = p
"""
    SimplePolynomial[{I}](p::AbstractPolynomialLike{C}, coefficient_type=C; kwargs...) where {C}

Creates a new `SimplePolynomial` based on any polynomial-like object that implements `MultivariatePolynomials`'s
`AbstractPolynomialLike` interface. The coefficients will be of type `coefficient_type`, the internal index type for the
monomials will be `I` (`UInt` if omitted). Keyword arguments are passed on to
[`SimpleMonomialVector`](@ref SimpleMonomialVector{I}(::AbstractVector{<:AbstractMonomialLike}, ::AbstractVector...)), which
allows to influence the variable mapping.
"""
function SimplePolynomial{I}(p::AbstractPolynomialLike{C1}, ::Type{C}=C1; kwargs...) where {C1,C,I<:Integer}
    coeffs = let c=coefficients(p)
        if eltype(c) == C
            collect(c)
        else
            convert(Vector{C}, collect(c))
        end
    end
    return SimplePolynomial(coeffs, SimpleMonomialVector{I}(monomials(p), coeffs; kwargs...))
end

SimplePolynomial(p::AbstractPolynomialLike, args...; kwargs...) = SimplePolynomial{UInt}(p, args...; kwargs...)

function (p::SimplePolynomial{C,Nr,Nc})(values::AbstractVector{V}) where {C,V,Nr,Nc}
    if length(values) == Nr + 2Nc
        all(v -> v[1] == conj(v[2]), zip(@view(values[Nr+1:Nr+Nc]), @view(values[Nr+Nc+1:end]))) ||
            throw(ArgumentError("Conjugate values did not match conjugate of their actual values"))
    else
        length(values) == Nr + Nc || throw(ArgumentError("Invalid number of values"))
    end
    T = promote_type(C, V)
    result = zero(T)
    real_values = @view(values[1:Nr])
    complex_values = @view(values[Nr+1:Nr+Nc])
    # This is not StaticPolynomials - we are interested in evaluating the polynomials, but only once or few times, so the
    # effort of generating efficient code is not really worth it.
    for t in p
        val = T(coefficient(t))
        for (var, pow) in monomial(t)
            idx = variable_index(var)
            if isreal(var)
                val *= real_values[idx]^pow
            elseif isconj(var)
                val *= conj(complex_values[idx])^pow
            else
                val *= complex_values[idx]^pow
            end
        end
        result += val
    end
    return result
end

function effective_variables_in(p::SimplePolynomial, in)
    for t in p
        iszero(coefficient(t)) || effective_variables_in(monomial(t), in) || return false
    end
    return true
end

"""
    change_backend(p::SimplePolynomial, variables::AbstractVector{<:AbstractVariable})

Changes a `SimplePolynomial` into a different implementation of `MultivariatePolynomials`, where the variables are taken from
the given vector in the order as they appear (but keeping real and complex variables distinct). Note that the coefficients are
returned without making a copy, which depending on the backend can imply back-effects on `p` itself should the output be
changed.
This conversion is not particularly efficient, as it works with generic implementations.
"""
change_backend(p::SimplePolynomial, variables::AbstractVector{<:AbstractVariable}) =
    polynomial(p.coeffs, change_backend(p.monomials, variables))
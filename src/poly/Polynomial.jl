export IntPolynomial

struct IntPolynomial{C,Nr,Nc,M<:IntMonomialVector{Nr,Nc}} <: AbstractPolynomial{C}
    coeffs::Vector{C}
    monomials::M

    """
        IntPolynomial(coeffs::Vector{C}, monomials::IntMonomialVector)

    Creates a `IntPolynomial` with well-defined coefficients and monomials, which must have the same length.
    """
    @inline function IntPolynomial(coeffs::Vector{C}, monomials::M) where {C,Nr,Nc,M<:IntMonomialVector{Nr,Nc}}
        length(coeffs) == length(monomials) || throw(ArgumentError("Lengths are different"))
        return new{C,Nr,Nc,M}(coeffs, monomials)
    end
end

MultivariatePolynomials.variables(p::XorTX{IntPolynomial}) = variables(monomial_type(p))

MultivariatePolynomials.coefficients(p::IntPolynomial) = p.coeffs

function MultivariatePolynomials.coefficient(p::IntPolynomial{C,Nr,Nc}, m::IntMonomial{Nr,Nc}) where {C,Nr,Nc}
    pos = searchsortedlast(p.monomials, m)
    (iszero(pos) || @inbounds(p.monomials[pos] != m)) && return zero(C)
    @inbounds return p.coeffs[pos]
end

MultivariatePolynomials.monomials(p::IntPolynomial) = p.monomials

MultivariatePolynomials.terms(p::IntPolynomial) = collect(p)

MultivariatePolynomials.nterms(p::IntPolynomial) = length(p.coeffs)

MultivariatePolynomials.nvariables(::XorTX{IntPolynomial{C,Nr,Nc} where {C}}) where {Nr,Nc} = Nr + 2Nc

Base.iterate(p::IntPolynomial) =
    isempty(p.coeffs) ? nothing : (@inbounds(term_type(p)(p.coeffs[begin], p.monomials[begin])), 2)
Base.iterate(p::IntPolynomial, i::Integer) =
    i > length(p.coeffs) ? nothing : (@inbounds(term_type(p)(p.coeffs[i], p.monomials[i])), i +1)
Base.IteratorSize(::Type{<:IntPolynomial}) = Base.HasLength()
Base.IteratorEltype(::Type{<:IntPolynomial}) = Base.HasEltype()
Base.eltype(::Type{P}) where {P<:IntPolynomial} = term_type(P)
Base.length(p::IntPolynomial) = length(p.coeffs)

for f in (:mindegree, :maxdegree, :extdegree, :mindegree_complex, :maxdegree_complex, :extdegree_complex, :minhalfdegree,
          :maxhalfdegree, :exthalfdegree, :variables)
    @eval MultivariatePolynomials.$f(p::IntPolynomial) = $f(p.monomials)
end

MultivariatePolynomials.constant_monomial(m::IntPolynomial) = constant_monomial(first(m.monomials))

MultivariatePolynomials.constant_monomial(P::Type{<:IntPolynomial}) = constant_monomial(eltype(P))

Base.zero(::Type{<:IntPolynomial{C,Nr,Nc,M}}) where {C,Nr,Nc,I<:Integer,E,M<:IntMonomialVector{Nr,Nc,I,E}} =
    IntPolynomial(
        C[],
        IntMonomialVector{Nr,Nc}(
            unsafe,
            ExponentsAll{Nr+2Nc,I}(),
            E <: Union{<:AbstractExponents,<:Tuple{AbstractExponents,AbstractUnitRange}} ? (one(I):zero(I)) : I[]
        )
    )

Base.zero(p::IntPolynomial{C,Nr,Nc,M}) where {C,Nr,Nc,I<:Integer,E,M<:IntMonomialVector{Nr,Nc,I,E}} =
    IntPolynomial(
        C[],
        IntMonomialVector{Nr,Nc}(
            unsafe,
            p.monomials.e,
            E <: Union{<:AbstractExponents,<:Tuple{AbstractExponents,AbstractUnitRange}} ? (one(I):zero(I)) : I[]
        )
    )

Base.one(::Type{<:IntPolynomial{C,Nr,Nc,M}}) where {C,Nr,Nc,I<:Integer,E,M<:IntMonomialVector{Nr,Nc,I,E}} =
    IntPolynomial(
        [one(C)],
        IntMonomialVector{Nr,Nc}(
            ExponentsAll{Nr+2Nc,I}(),
            E <: Union{<:AbstractExponents,<:Tuple{AbstractExponents,AbstractUnitRange}} ? (one(I):one(I)) : [one(I)]
        )
    )

function Base.one(p::IntPolynomial{C,Nr,Nc,M}) where {C,Nr,Nc,I<:Integer,E,M<:IntMonomialVector{Nr,Nc,I,E}}
    eind = convert_index(p.monomials.e, ExponentsAll{Nr+2Nc,I}(), one(I))
    return IntPolynomial(
        [one(C)],
        IntMonomialVector{Nr,Nc}(
            unsafe,
            p.monomials.e,
            E <: Union{<:AbstractExponents,<:Tuple{AbstractExponents,AbstractUnitRange}} ? (eind:eind) : [eind]
        )
    )
end

MultivariatePolynomials.map_coefficients(f::Function, p::IntPolynomial) = IntPolynomial(map(f, p.coeffs), p.monomials)

function MultivariatePolynomials.map_coefficients!(f::Function, p::IntPolynomial)
    @inbounds for i in 1:length(p.coeffs)
        p.coeffs[i] = f(p.coeffs[i])
    end
    return p
end

function MultivariatePolynomials.effective_variables(p::IntPolynomial{<:Any,Nr,Nc}; rettype::Type{V}=Vector, by::F=identity) where {Nr,Nc,V<:Union{Vector,Set},F}
    ev = Set{Base.promote_op(by, IntVariable{Nr,Nc,smallest_unsigned(Nr+2Nc)})}()
    for t in p, (v, _) in monomial(t)
        push!(ev, by(v))
    end
    V <: Vector && return sort!(collect(ev))
    return ev
end

# these are tricky. They shouldn't even be there since they might (or will) violate the no-allocation policy of
# IntPolynomials
function Base.conj(p::IntPolynomial)
    newcoeffs = conj(p.coeffs)
    IntPolynomial(newcoeffs, conj(p.monomials, newcoeffs)) # conjugation implies reordering, so sort the coefficients along
end
function MultivariatePolynomials.LinearAlgebra.adjoint(p::IntPolynomial)
    newcoeffs = adjoint.(p.coeffs)
    IntPolynomial(newcoeffs, conj(p.monomials, newcoeffs))
end
Base.isreal(p::IntPolynomial{<:Any,<:Any,0}) = all(∘(iszero, imag), p.coeffs)
function Base.isreal(p::IntPolynomial)
    # the monomial vector is sorted and unique
    for t in p
        coefficient(t) == conj(coefficient(p, IntConjMonomial(monomial(t)))) || return false
    end
    return true
end

MultivariatePolynomials.LinearAlgebra.transpose(p::IntPolynomial) = IntPolynomial(transpose.(p.coeffs), p.monomials)

IntPolynomial(p::IntPolynomial) = p
"""
    IntPolynomial[{I}](p::AbstractPolynomialLike{C}, coefficient_type=C; kwargs...) where {C}

Creates a new `IntPolynomial` based on any polynomial-like object that implements `MultivariatePolynomials`'s
`AbstractPolynomialLike` interface. The coefficients will be of type `coefficient_type`, the internal index type for the
monomials will be `I` (`UInt` if omitted). Keyword arguments are passed on to
[`IntMonomialVector`](@ref IntMonomialVector{I}(::AbstractVector{<:AbstractMonomialLike}, ::AbstractVector...)), which
allows to influence the variable mapping.
"""
function IntPolynomial{I}(p::AbstractPolynomialLike{C1}, ::Type{C}=C1; kwargs...) where {C1,C,I<:Integer}
    coeffs = let c=coefficients(p)
        if eltype(c) == C
            collect(c)
        else
            convert(Vector{C}, collect(c))
        end
    end
    return IntPolynomial(coeffs, IntMonomialVector{I}(monomials(p), coeffs; kwargs...))
end

IntPolynomial(p::AbstractPolynomialLike, args...; kwargs...) = IntPolynomial{UInt}(p, args...; kwargs...)

function (p::IntPolynomial{C,Nr,Nc})(values::AbstractVector{V}) where {C,V,Nr,Nc}
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

function effective_variables_in(p::IntPolynomial, in)
    for t in p
        iszero(coefficient(t)) || effective_variables_in(monomial(t), in) || return false
    end
    return true
end

function combine_polynomials(p₁::IntPolynomial{C1,Nr,Nc}, p₂::IntPolynomial{C2,Nr,Nc},
    combine::Union{typeof(+),typeof(-)}) where {C1,C2,Nr,Nc}
    # This is in principle similar to merge_monomial_vectors, but we can do a much simpler implementation. We will assume that
    # the exponent sets are the same
    m₁, m₂ = monomials(p₁), monomials(p₂)
    E₁, E₂ = m₁.e, m₂.e
    E₁ == E₂ || throw(ArgumentError("Combination of polynomials is implemented only for the same exponent sets."))
    if isempty(m₁)
        return IntPolynomial(x -> combine(zero(C1), x), coefficients(p₂), m₂)
    elseif isempty(m₂)
        return IntPolynomial(x -> combine(x, zero(C2)), coefficients(p₁), m₁)
    elseif m₁ == m₂
        return IntPolynomial(map(splat(combine), zip(coefficients(p₁), coefficients(p₂))), m₁)
    elseif last(m₁.indices) < first(m₂.indices)
        newcoeffs = Vector{promote_type(C1, C2)}(undef, length(m₁) + length(m₂))
        map!(x -> combine(x, zero(C2)), @view(newcoeffs[1:length(m₁)]), coefficients(p₁))
        map!(x -> combine(zero(C1), x), @view(newcoeffs[length(m₁)+1:end]), coefficients(p₂))
        if m₁.indices isa AbstractUnitRange && m₂.indices isa AbstractUnitRange && last(m₁.indices) == first(m₂.indices) -1
            return IntPolynomial(newcoeffs, first(m₁.indices):last(m₂.indices))
        else
            return IntPolynomial(newcoeffs, [m₁.indices; m₂.indices])
        end
    elseif last(m₂.indices) < first(m₁.indices)
        newcoeffs = Vector{promote_type(C1, C2)}(undef, length(m₁) + length(m₂))
        map!(x -> combine(zero(C1), x), @view(newcoeffs[1:length(m₂)]), coefficients(p₂))
        map!(x -> combine(x, zero(C2)), @view(newcoeffs[length(m₂)+1:end]), coefficients(p₁))
        if m₁.indices isa AbstractUnitRange && m₂.indices isa AbstractUnitRange && last(m₂.indices) == first(m₁.indices) -1
            return IntPolynomial(newcoeffs, first(m₂.indices):last(m₁.indices))
        else
            return IntPolynomial(newcoeffs, [m₁.indices; m₂.indices])
        end
    elseif m₁.indices isa AbstractUnitRange && m₂.indices isa AbstractUnitRange &&
        ((first(m₁.indices) ≤ last(m₂.indices) && last(m₁.indices) ≥ first(m₂.indices)) ||
            (first(m₂.indices) ≤ last(m₂.indices) && last(m₂.indices) ≤ first(m₁.indices)))
        newcoeffs = Vector{promote_type(C1, C2)}(undef, length(m₁) + length(m₂))
        indices₁, indices₂ = m₁.indices, m₂.indices
        coeffs₁, coeffs₂ = coefficients(p₁), coefficients(p₂)
        resize!(newcoeffs, Solver.count_uniques(indices₁, indices₂, (i, i₁, i₂) -> @inbounds begin
            if ismissing(i₁)
                newcoeffs[i] = combine(zero(C1), coeffs₂[i₂])
            else
                newcoeffs[i] = combine(coeffs₁[i₁], ismissing(i₂) ? zero(C2) : coeffs₂[i₂])
            end
        end))
        return IntPolynomial(newcoeffs, min(first(m₁.indices), first(m₂.indices)):max(last(m₁.indices), last(m₂.indices)))
    else
        newcoeffs = Vector{promote_type(C1, C2)}(undef, length(m₁) + length(m₂))
        newindices = similar(newcoeffs, _get_I(E₁))
        indices₁, indices₂ = m₁.indices, m₂.indices
        coeffs₁, coeffs₂ = coefficients(p₁), coefficients(p₂)
        resize!(newcoeffs, Solver.count_uniques(indices₁, indices₂, (i, i₁, i₂) -> @inbounds begin
            if ismissing(i₁)
                newcoeffs[i] = combine(zero(C1), coeffs₂[i₂])
                newindices[i] = indices₂[i₂]
            else
                newcoeffs[i] = combine(coeffs₁[i₁], ismissing(i₂) ? zero(C2) : coeffs₂[i₂])
                newindices[i] = indices₁[i₁]
            end
        end))
        resize!(newindices, length(newcoeffs))
        return IntPolynomial(newcoeffs, newindices)
    end
end

"""
    change_backend(p::IntPolynomial, variables::AbstractVector{<:AbstractVariable})

Changes a `IntPolynomial` into a different implementation of `MultivariatePolynomials`, where the variables are taken from
the given vector in the order as they appear (but keeping real and complex variables distinct). Note that the coefficients are
returned without making a copy, which depending on the backend can imply back-effects on `p` itself should the output be
changed.
This conversion is not particularly efficient, as it works with generic implementations.
"""
change_backend(p::IntPolynomial, variables::AbstractVector{<:AbstractVariable}) =
    polynomial(p.coeffs, change_backend(p.monomials, variables))
# This file is only loaded if Oscar is present and adapts all necessary functions to also provide implementations for Oscar
# polynomials

import .Oscar

# conveniences for identifying a ring with its list of generators
Base.length(r::Oscar.Ring) = Oscar.ngens(r)
# ordinary MP functions. Note that in Oscar, we do not distinguish variables, monomials, terms, and polynomials
function MultivariatePolynomials.variable(p::Oscar.MPolyElem{T}) where {T}
    Oscar.is_gen(p) || throw(InexactError(:variable, AbstractVariable, p))
    return p
end
function MultivariatePolynomials.name(p::Oscar.MPolyElem{T}) where {T}
    io = IOBuffer()
    print(io, variable(p))
    seekstart(io)
    return read(io, String)
end
MultivariatePolynomials.variable_union_type(::Union{T,Type{T}}) where {Q,T<:Oscar.MPolyElem{Q}} = T
function MultivariatePolynomials.monomial_type(m::T) where {Q,T<:Oscar.MPolyElem{Q}}
    @assert(Oscar.is_monomial(m))
    return T
end
MultivariatePolynomials.monomial_type(::Type{T}) where {Q,T<:Oscar.MPolyElem{Q}} = T
MultivariatePolynomials.variables(p::Oscar.MPolyElem{T}) where {T} = Oscar.gens(Oscar.parent(p))
MultivariatePolynomials.effective_variables(p::Oscar.MPolyElem{T}) where {T} = Oscar.vars(p)
MultivariatePolynomials.nvariables(p::Oscar.MPolyElem{T}) where {T} = Oscar.nvars(Oscar.parent(p))
function MultivariatePolynomials.exponents(p::Oscar.MPolyElem{T}) where {T}
    @assert(Oscar.is_term(p))
    return Oscar.exponent_vector(p, 1)
end
function MultivariatePolynomials.degree(p::Oscar.MPolyElem{T}) where {T}
    @assert(Oscar.is_term(p))
    return Oscar.total_degree(p)
end
function MultivariatePolynomials.degree(p::T, v::T) where {Q,T<:Oscar.MPolyElem{Q}}
    @assert(Oscar.is_term(p) && Oscar.parent(p) === Oscar.parent(v))
    return Oscar.degree(p, v)
end
function MultivariatePolynomials.isconstant(p::Oscar.MPolyElem{T}) where {T}
    @assert(Oscar.is_term(p))
    return Oscar.is_constant(p)
end
MultivariatePolynomials.powers(p::Oscar.MPolyElem{T}) where {T} =
    ((v, e) for (v, e) in zip(variables(p), exponents(p)) if !iszero(e))
MultivariatePolynomials.constant_monomial(p::Oscar.MPolyElem{T}) where {T} = one(p)
function MultivariatePolynomials.map_exponents(f, m1::M, m2::M) where {Q,M<:Oscar.MPolyElem{Q}}
    ring = Oscar.parent(m1)
    @assert(Oscar.is_monomial(m1) && Oscar.is_monomial(m2) && Oscar.parent(m2) === ring)
    return ring([one(Q)], Vector{Int}[[f(α, β) for (α, β) in zip(Oscar.exponent_vector(m1, 1), Oscar.exponent_vector(m2, 1))]])
end
function MultivariatePolynomials.term(coef::T, p::Oscar.MPolyElem{T}) where {T}
    # we always copy, else we would have to mutate p
    @assert(Oscar.is_monomial(p))
    return coef * p
end
MultivariatePolynomials.term_type(::Union{T,Type{T}}) where {Q,T<:Oscar.MPolyElem{Q}} = T
# MultivariatePolynomials.term_type(::Union{T,Type{T}}, ::Type{C}) where {Q,T<:Oscar.MPolyElem{Q},C} --- would need to be more specific, not just MPolyElem but probably Nemo.fmpq_mpoly, breaking genericity
function MultivariatePolynomials.coefficient(t::Oscar.MPolyElem{T}) where {T}
    @assert(Oscar.is_term(t))
    return Oscar.coeff(t, 1)
end
function MultivariatePolynomials.coefficient(p::T, m::T) where {Q,T<:Oscar.MPolyElem{Q}}
    @assert(Oscar.is_monomial(m) && Oscar.parent(p) === Oscar.parent(m))
    return Oscar.coeff(p, Oscar.exponent_vector(m, 1))
end
MultivariatePolynomials.coefficient_type(::Union{T,Type{T}}) where {Q,T<:Oscar.MPolyElem{Q}} = Q
function MultivariatePolynomials.monomial(t::Oscar.MPolyElem{T}) where {T}
    @assert(Oscar.is_term(t))
    return Oscar.monomial(t, 1)
end
# we cannot implement the following two for the type version
MultivariatePolynomials.constant_term(α::T, p::Oscar.MPolyElem{T}) where {T} = α * one(p)
MultivariatePolynomials.zero_term(p::Oscar.MPolyElem{T}) where {T} = zero(p)
MultivariatePolynomials.polynomial(p::Oscar.MPolyElem{T}) where {T} = p
# MultivariatePolynomials.polynomial(p::Oscar.MPolyElem{T}, ::Type{C}) where {T,C} --- would need to create new ring
function MultivariatePolynomials.polynomial(a::AbstractVector{Q}, mv::AbstractVector{T}) where {Q,T<:Oscar.MPolyElem{Q}}
    @assert(length(a) == length(mv))
    isempty(mv) && return zero(T)
    ring = Oscar.parent(first(mv))
    @assert(all(Oscar.is_monomial.(mv)) && all(Oscar.parent.(mv) .=== (ring,)))
    return ring(collect(a), Oscar.exponent_vector.(mv, 1))
end
MultivariatePolynomials.polynomial(terms::AbstractVector{T},
    ::MultivariatePolynomials.ListState=MultivariatePolynomials.MessyState()) where {Q,T<:Oscar.MPolyElem{Q}} = sum(terms)
function MultivariatePolynomials.polynomial(f::Function, mv::AbstractVector{T}) where {Q,T<:Oscar.MPolyElem{Q}}
    isempty(mv) && return zero(T)
    ring = Oscar.parent(first(mv))
    @assert(all(Oscar.is_monomial.(mv)) && all(Oscar.parent.(mv) .=== (ring,)))
    return ring(f.(1:length(mv)), Oscar.exponent_vector.(mv, 1))
end
MultivariatePolynomials.polynomial_type(::Union{T,Type{T}}) where {Q,T<:Oscar.MPolyElem{Q}} = T
MultivariatePolynomials.polynomial_type(::Union{T,Type{T}}, ::Type{C}) where {Q,T<:Oscar.MPolyElem{Q},C} = Oscar.MPolyElem{C}
MultivariatePolynomials.terms(p::Oscar.MPolyElem{T}) where {T} = Oscar.terms(p)
Base.iterate(p::Oscar.MPolyElem{T}, args...) where {T} = iterate(Oscar.terms(p), args...)
MultivariatePolynomials.nterms(p::Oscar.MPolyElem{T}) where {T} = length(p)
MultivariatePolynomials.coefficients(p::Oscar.MPolyElem{T}) where {T} = Oscar.coefficients(p)
MultivariatePolynomials.coefficients(p::T, X::AbstractVector{T}) where {Q,T<:Oscar.MPolyElem{Q}} =
    (Oscar.coeff(p, m) for m in X)
# MultivariatePolynomials.coefficient(p::T, m::T, vars) --- not compatible with Oscar's point of view
MultivariatePolynomials.monomials(p::Oscar.MPolyElem{T}) where {T} = Oscar.monomials(p)
function fillZfordeg!(Z, ring::R, n, deg, filter::Function) where {R<:Oscar.AbstractAlgebra.Ring}
    # from  DynamicPolynomials
    z = zeros(Int, n)
    z[1] = deg
    while true
        mon = ring([one(coefficient_type(Oscar.elem_type(ring)))], [z])
        filter(mon) && push!(Z, mon)
        z[end] == deg && break
        sum = 1
        for j in (n-1):-1:1
            if z[j] != 0
                z[j] -= 1
                z[j+1] += sum
                break
            else
                sum += z[j+1]
                z[j+1] = 0
            end
        end
    end
end
function MultivariatePolynomials.monomials(ring::Oscar.Ring, degs::AbstractVector{Int}, filter::Function=m->true)
    Z = Vector{typeof(zero(ring))}()
    n = Oscar.nvars(ring)
    for deg in sort(degs, rev=true)
        fillZfordeg!(Z, ring, n, deg, filter)
    end
    return Z
end
MultivariatePolynomials.ordering(p::Oscar.MPolyElem{T}) where {T} = Oscar.ordering(Oscar.parent(p))
MultivariatePolynomials.mindegree(X::AbstractVector{<:Oscar.MPolyElem{T}}) where {T} =
    isempty(X) ? 0 : minimum(t -> Oscar.total_degree(t), X)
MultivariatePolynomials.mindegree(X::AbstractVector{<:T}, v::T) where {Q,T<:Oscar.MPolyElem{Q}} =
    isempty(X) ? 0 : minimum(t -> Oscar.degree(t, v), X)
MultivariatePolynomials.mindegree(p::Oscar.MPolyElem{T}, args...) where {T} = mindegree(collect(Oscar.monomials(p)), args...)
MultivariatePolynomials.maxdegree(X::AbstractVector{<:Oscar.MPolyElem{T}}) where {T} =
    mapreduce(t -> Oscar.total_degree(t), max, X, init=0)
MultivariatePolynomials.maxdegree(X::AbstractVector{<:Oscar.MPolyElem{T}}, v::Oscar.MPolyElem{T}) where {T} =
    mapreduce(t -> Oscar.degree(t, v), max, X, init=0)
MultivariatePolynomials.maxdegree(p::Oscar.MPolyElem{T}, args...) where {T} = maxdegree(collect(Oscar.monomials(p)), args...)
MultivariatePolynomials.extdegree(p::Union{T,AbstractVector{<:T}}, args...) where {Q,T<:Oscar.MPolyElem{Q}} = (mindegree(p, args...), maxdegree(p, args...))
MultivariatePolynomials.leading_term(p::Oscar.MPolyElem{T}) where {T} = Oscar.leading_term(p)
MultivariatePolynomials.leading_coefficient(p::Oscar.MPolyElem{T}) where {T} = Oscar.leading_coefficient(p)
MultivariatePolynomials.leading_monomial(p::Oscar.MPolyElem{T}) where {T} = Oscar.leading_monomial(p)
MultivariatePolynomials.remove_leading_term(p::Oscar.MPolyElem{T}) where {T} = Oscar.tail(p)
function MultivariatePolynomials.remove_monomials(p::T, mons::AbstractArray{T}) where {Q,T<:Oscar.MPolyElem{Q}}
    isempty(mons) && return p
    ring = Oscar.parent(p)
    @assert(all(Oscar.parent.(mons) .=== (ring,)))
    result = copy(p)
    Oscar.setcoeff!.(result, Oscar.exponent_vector.(mons, 1), zero(Q))
    return result
end
MultivariatePolynomials.monic(p::Oscar.MPolyElem{T}) where {T} = Oscar.divexact(p, Oscar.leading_coefficient(p))
MultivariatePolynomials.map_coefficients(f::Function, p::Oscar.MPolyElem{T}, ::Bool=false) where {T} = Oscar.map_coefficients(f, p)
function MultivariatePolynomials.map_coefficients!(f::Function, p::Oscar.MPolyElem{T}, ::Bool=false) where {T}
    for i in 1:length(p)
        Oscar.setcoeff!(p, i, f(Oscar.coeff(p, i)))
    end
    return p
end
function MultivariatePolynomials.map_coefficients_to!(output::T, f::Function, p::T, ::Bool=false) where {Q,T<:Oscar.MPolyElem{Q}}
    for i in length(output):-1:1
        Oscar.setcoeff!(p, i, zero(Q))
    end
    for i in 1:length(p)
        Oscar.setcoeff!(p, Oscar.exponent_vector(p, i), f(Oscar.coeff(p, i)))
    end
    return output
end
function MultivariatePolynomials.monomial_vector(X::AbstractVector{MT}) where {T,MT<:Oscar.MPolyElem{T}}
    isempty(X) && return X
    ring = Oscar.parent(first(X))
    @assert(all(Oscar.is_monomial.(X)) && all(Oscar.parent.(X) .=== (ring,)))
    Y = sort(X, rev=true)
    dups = findall(i -> Y[i] == Y[i-1], 2:length(Y))
    deleteat!(Y, dups)
    return Y
end
function MultivariatePolynomials.monomial_vector(a::AbstractVector, X::AbstractVector{MT}) where {T,MT<:Oscar.MPolyElem{T}}
    @assert(length(a) == length(X))
    isempty(X) && return X
    ring = Oscar.parent(first(X))
    @assert(all(Oscar.is_monomial.(X)) && all(Oscar.parent.(X) .=== (ring,)))
    σ = sortperm(X, rev=true)
    b = T[]
    sizehint!(b, length(a))
    Y = MT[]
    sizehint!(Y, length(X))
    push!(b, zero(T))
    push!(Y, X[first(σ)])
    for (coeff, x) in zip(@view(a[σ]), @view(X[σ]))
        if x == Y[end]
            b[end] += coeff
        else
            push!(b, coeff)
            push!(Y, x)
        end
    end
    sizehint!(b, length(b))
    sizehint!(Y, length(Y))
    return b, Y
end
MultivariatePolynomials.monomial_vector_type(::AbstractVector{MT}) where {Q,MT<:Oscar.MPolyElem{Q}} = Vector{MT}
MultivariatePolynomials.empty_monomial_vector(::T) where {Q,T<:Oscar.MPolyElem{Q}} = T[]
function MultivariatePolynomials.sort_monomial_vector(X::AbstractVector{MT}) where {Q,MT<:Oscar.MPolyElem{Q}}
    isempty(X) && return Int[], MT[]
    ring = Oscar.parent(first(X))
    @assert(all(Oscar.is_monomial.(X)) && all(Oscar.parent.(X) .=== (ring,)))
    σ = sortperm(X, rev=true)
    dups = findall(i -> X[σ[i]] == X[σ[i-1]], 2:length(σ))
    deleteat!(σ, dups)
    return σ, X[σ]
end
MultivariatePolynomials.merge_monomial_vectors(ms::AbstractVector{MVT}) where {Q,MT<:Oscar.MPolyElem{Q},MVT<:AbstractVector{MT}} =
    sort_monomial_vector(vcat(ms...))[2]
# ordinary division and multiplication
Base.:*(M::Oscar.MPolyElem{T}, v::AbstractFloat) where {T} = M * rationalize(v)
Base.:*(v::AbstractFloat, M::Oscar.MPolyElem{T}) where {T} = rationalize(v) * M
Base.:/(M::Oscar.MPolyElem{T}, v::Union{Integer,Rational}) where {T} = M // v
Base.:/(M::Oscar.MPolyElem{T}, v::AbstractFloat) where {T} = M // rationalize(v)
# substitution
const OscarSubstitution{T,V<:Oscar.MPolyElem{T}} = Pair{V}
const OscarMultiSubstitution{T,V<:Oscar.MPolyElem{T},N} = Pair{<:NTuple{N,V}, <:NTuple{N,Any}}
const OscarMultiVectorSubstitution{T,V<:Oscar.MPolyElem{T}} = Pair{<:Tuple{Vararg{V}}, <:AbstractVector}
const OscarVectorMultiSubstitution{T,V<:Oscar.MPolyElem{T}} = Pair{<:AbstractVector{V}, <:Tuple}
const OscarVectorMultiVectorSubstitution{T,V<:Oscar.MPolyElem{T}} = Pair{<:AbstractVector{V}, <:AbstractVector}
const OscarAbstractMultiSubstitution{T,V<:Oscar.MPolyElem{T}} = Union{OscarMultiSubstitution{T,V}, OscarMultiVectorSubstitution{T,V}, OscarVectorMultiSubstitution{T,V}, OscarVectorMultiVectorSubstitution{T,V}}
const OscarAbstractSubstitution{T,V<:Oscar.MPolyElem{T}} = Union{OscarSubstitution{T,V}, OscarAbstractMultiSubstitution{T,V}}

MultivariatePolynomials.subs(p::M, s::OscarAbstractSubstitution{T,M}...) where {T,M<:Oscar.MPolyElem{T}} = subs(p, M[], [], s...)
function MultivariatePolynomials.subs(p::M, vars::Vector{M}, vals::Vector) where {T,M<:Oscar.MPolyElem{T}}
    typed_vals = identity.(vals)
    if eltype(typed_vals) <: AbstractFloat
        return eltype(typed_vals).(Oscar.evaluate(p, vars, rationalize.(typed_vals)))
    else
        return Oscar.evaluate(p, vars, typed_vals)
    end
end
function MultivariatePolynomials.subs(p::M, vars::Vector{M}, vals::Vector, substitution::OscarSubstitution{T,M}, s::OscarAbstractSubstitution{T,M}...) where {T,M<:Oscar.MPolyElem{T}}
    @assert(Oscar.is_gen(substitution.first) && Oscar.parent(p) == Oscar.parent(substitution.first))
    push!(vars, substitution.first)
    push!(vals, substitution.second)
    return subs(p, vars, vals, s...)
end
function MultivariatePolynomials.subs(p::M, vars::Vector{M}, vals::Vector, substitution::OscarAbstractMultiSubstitution{T,M}, s::OscarAbstractSubstitution{T,M}...) where {T,M<:Oscar.MPolyElem{T}}
    @assert(length(substitution.first) == length(substitution.second) && all(Oscar.is_gen.(substitution.first)) &&
        all(Oscar.parent.(substitution.first) .== Oscar.parent(p)))
    append!(vars, substitution.first)
    append!(vals, substitution.second)
    return subs(p, vars, vars, s...)
end
# evaluation is already defined in Oscar, but only in a very specific way. MultivariatePolynomials is more general
function (x::Oscar.MPolyElem{T})(s::OscarAbstractMultiSubstitution{T,M}) where {T,M<:Oscar.MPolyElem{T}}
    @assert(all(Oscar.is_gen.(s.first)))
    length(s.first) == Oscar.nvars(Oscar.parent(x)) || throw(ArgumentError("Not all variables were assigned a value. Use `subs` to substitute only a subset of the variables."))
    σ = sortperm(s.first, rev=true)
    print(σ)
    isnothing(findfirst(i -> s.first[σ[i]] == s.first[σ[i-1]], 2:length(σ))) || throw(ArgumentError("Not all variables were assigned a value. Use `subs` to substitute only a subset of the variables."))
    return x(@view(s.second[σ])...)
end
# evaluation is already defined in Oscar, but only in a very specific way. MultivariatePolynomials is more general
function (x::Oscar.MPolyElem{T})(s::OscarAbstractMultiSubstitution{T,M}) where {T,M<:Oscar.MPolyElem{T}}
    @assert(typeof(x) === M)
    length(s.first) == Oscar.nvars(Oscar.parent(x)) || throw(ArgumentError("Not all variables were assigned a value. Use `subs` to substitute only a subset of the variables."))
    σ = sortperm(s.first, rev=true)
    isnothing(findfirst(i -> s.first[σ[i]] == s.first[σ[i-1]], 2:length(σ))) || throw(ArgumentError("Not all variables were assigned a value. Use `subs` to substitute only a subset of the variables."))
    return x(@view(s.second[σ])...)
end
# division
MultivariatePolynomials.divides(t1::T, t2::T) where {Q,T<:Oscar.MPolyElem{Q}} = Oscar.divides(t2, t1)[1]
Base.rem(p::T, gb::Oscar.MPolyIdeal{T}) where {Q,T<:Oscar.MPolyElem{Q}} = Oscar.normal_form(p, gb)
Base.rem(p::T, ::EmptyGröbnerBasis) where {Q,T<:Oscar.MPolyElem{Q}} = p

# matrix polynomials
MultivariatePolynomials.variables(m::AbstractMatrix{T}) where {Q,T<:Oscar.MPolyElem{Q}} = union(variables.(m)...)

MultivariatePolynomials.monomials(m::AbstractMatrix{T}) where {Q,T<:Oscar.MPolyElem{Q}} = union(monomials.(m))

MultivariatePolynomials.coefficients(m::AbstractMatrix{T}) where {Q,T<:Oscar.MPolyElem{Q}} = coefficients.(m, (monomials(m),))
MultivariatePolynomials.coefficients(m::AbstractMatrix{T}, X::AbstractVector) where {Q,T<:Oscar.MPolyElem{Q}} = coefficients.(m, (X,))

MultivariatePolynomials.mindegree(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}} = minimum(mindegree(p, args...) for p in m)
mindegree_complex(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}} = minimum(mindegree_complex(p, args...) for p in m)
minhalfdegree(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}} = minimum(minhalfdegree(p, args...) for p in m)
MultivariatePolynomials.maxdegree(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}} = maximum(maxdegree(p, args...) for p in m)
maxdegree_complex(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}} = maximum(maxdegree_complex(p, args...) for p in m)
maxhalfdegree(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}} = maximum(maxhalfdegree(p, args...) for p in m)
function MultivariatePolynomials.extdegree(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}}
    l = typemax(Int)
    u = 0
    for p in m
        (newl, newu) = extdegree(p, args...)
        newl < l && (l = newl)
        newu > u && (u = newu)
    end
    return l, u
end
function extdegree_complex(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}}
    l = typemax(Int)
    u = 0
    for p in m
        (newl, newu) = extdegree_complex(p, args...)
        newl < l && (l = newl)
        newu > u && (u = newu)
    end
    return l, u
end
function exthalfdegree(m::AbstractMatrix{T}, args...) where {Q,T<:Oscar.MPolyElem{Q}}
    l = typemax(Int)
    u = 0
    for p in m
        (newl, newu) = exthalfdegree(p, args...)
        newl < l && (l = newl)
        newu > u && (u = newu)
    end
    return l, u
end

function MultivariatePolynomials.divides(t1s::AbstractVector{T}, t2) where {Q,T<:Oscar.MPolyElem{Q}}
    for t1 in t1s
        if divides(t1, t2)
            return true
        end
    end
    return false
end

MultivariatePolynomials.effective_variables(m::AbstractMatrix{T}) where {Q,T<:Oscar.MPolyElem{Q}} = union(effective_variables.(m)...)

# complex MP functions. We do not implement this, only define the fallbacks.
Base.isreal(::Oscar.MPolyElem{T}) where {T} = true
isrealpart(::Oscar.MPolyElem{T}) where {T} = false
isimagpart(::Oscar.MPolyElem{T}) where {T} = false
isconj(::Oscar.MPolyElem{T}) where {T} = false
ordinary_variable(x::Oscar.MPolyElem{T}) where {T} = x
ordinary_variable(x::AbstractVector{<:Oscar.MPolyElem{T}}) where {T} = copy(x)
Base.conj(x::Oscar.MPolyElem{T}) where {T} = x
Base.real(x::Oscar.MPolyElem{T}) where {T} = x
Base.imag(x::Oscar.MPolyElem{T}) where {T} = zero(x)
LinearAlgebra.adjoint(x::Oscar.MPolyElem{T}) where {T} = x
degree_complex(p::Oscar.MPolyElem{T}, args...) where {T} = degree(p, args...)
halfdegree(p::Oscar.MPolyElem{T}) where {T} = degree(p)
mindegree_complex(p::Oscar.MPolyElem{T}, args...) where {T} = mindegree(p, args...)
minhalfdegree(p::Oscar.MPolyElem{T}, args...) where {T} = mindegree(p, args...) ÷ 2
maxdegree_complex(p::Oscar.MPolyElem{T}, args...) where {T} = maxdegree(p, args...)
maxhalfdegree(p::Oscar.MPolyElem{T}, args...) where {T} = maxdegree(p, args...) ÷ 2
# we don't define MonomialComplexContainer, as Oscar polynomials cannot be complex

function poly_problem(objective::P, degree::Int; zero::AbstractVector{P}=P[],
    nonneg::AbstractVector{P}=P[], psd::AbstractVector{<:AbstractMatrix{P}}=AbstractMatrix{P}[],
    custom_basis::AbstractVector{P}=P[], perturbation::Union{Float64,<:AbstractVector{Float64}}=0.0) where {T,P<:Oscar.MPolyElem{T}}
    @assert all(Oscar.is_monomial.(custom_basis))
    ring = Oscar.parent(objective)
    @assert(all(Oscar.parent.(zero) .== (ring,)) && all(Oscar.parent.(nonneg) .== (ring,)) &&
        all(Oscar.parent.(psd) .== (ring,)) && all(Oscar.parent.(custom_basis) .== (ring,)))
    return poly_problem(objective, Oscar.gens(ring), degree, zero, nonneg, psd, custom_basis, perturbation)
end

function subbasis(b::Oscar.MPolyIdeal{T}, variables) where {T}
    @assert(!complex)
    gens = Oscar.gens(b)
    return Oscar.ideal(Oscar.base_ring(b), gens[effective_variables.(gens).⊆(variables,)])
end

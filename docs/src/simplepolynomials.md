```@meta
CurrentModule = PolynomialOptimization.SimplePolynomials.MultivariateExponents
```

# SimplePolynomials
`PolynomialOptimization` allows data in the form of any implementation that supports the
[`MultivariatePolynomials`](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) interface. However, it does not keep
the data in this way, which would not be particularly efficient. Instead, it is converted into an internal format, the
`SimplePolynomial`. It offers very compact storage (in some way even more compact than
[`SIMDPolynomials`](https://github.com/YingboMa/SIMDPolynomials.jl)) and is particularly focused on being as
allocation-free as possible - which means that once the original polynomials were created, at no further stage in processing
the polynomials or monomial bases will any allocations be done.

## AbstractExponents
For this, we first define a basic set of exponents which can occur in our monomial; this is a subtype of
[`AbstractExponents`](@ref). Such a set of exponents is always parameterized by the number of variables and by the integer data
type that contains the index of the monomials in the set. By default, this is `UInt`, but any Integer descendant can be chosen.
Note that with ``n`` variables (complex and conjugate variables counted separately) of maximum degree ``d``, the datatype must
be able to hold the value ``\binom{n + d}{n}``.
In the typical scenario with not too high degrees, machine data types should be sufficient (say, ``d = 4``, then even `UInt32`
can hold up to 564 variables, and with `UInt64`, more than ``10^5`` variables are possible).
A monomial vector is then either a complete cover of a degree-bound exponent set or a finite subset of any exponent set. It
needs no extra space apart from the description of the exponent set itself if it is a complete cover, and only the space
required to describe the subindices (they need not necessarily be a vector, they can also be a range) if it is a subset.

## Exponent set types
```@docs
AbstractExponents
AbstractExponentsUnbounded
AbstractExponentsDegreeBounded
ExponentsAll
ExponentsDegree
ExponentsMultideg
indextype
```

## Working with exponents
Exponent sets can be indexed and iterated. If consecutive elements are required, iteration is slightly faster, as indexing
requires to determine the degree for every operation. However, in this case it is fastest to preallocate a vector that can hold
the exponents and [`iterate!`](@ref iterate!(::AbstractVector{Int}, ::AbstractExponents)) with mutation of this vector, or use
the [`veciter`](@ref) wrapper for this task:
```@docs
iterate!(::AbstractVector{Int}, ::AbstractExponents)
veciter(::AbstractExponents, ::AbstractVector{Int})
```

When working with individual indices or exponents, conversion functions are provided.
```@docs
exponents_to_index
exponents_from_index(::AbstractExponents{<:Any,I}, ::I) where {I<:Integer}
exponents_sum
exponents_product
```
Note that `exponents_from_index` returns a lazy implementation of an `AbstractVector{Int}`; if the same exponents must be
accessed multple times, it might be beneficial to `collect` the result or copy it to a pre-allocated vector.

Further information can be obtained about one or two indices:
```@docs
degree_from_index(::AbstractExponents{N,I}, ::I) where {N,I<:Integer}
convert_index(::AbstractExponents{N,I}, ::AbstractExponents{N,IS}, ::IS, ::Int) where {N,I<:Integer,IS<:Integer}
compare_indices(::AbstractExponents{N,I1}, ::I1, ::_CompareOp, ::AbstractExponents{N,I2}, ::I2, ::Int) where {N,I1<:Integer,I2<:Integer}
```

Degree-bound exponent sets have a length:
```@docs
Base.length(::AbstractExponentsDegreeBounded)
```

## Internals
The internal cache of exponents is extremely important for fast access. Unless the unsafe versions of the functions are used,
it is always ensured that the cache is large enough to do the required operations; this is a quick check, but it can be elided
using the unsafe versions. They, as well as the functions that allow direct access to the cache, must only be used if you know
exactly what you are doing.
```@docs
index_counts
exponents_from_index(::Unsafe, ::AbstractExponents{<:Any,I}, ::I, ::Int) where {I<:Integer}
degree_from_index(::Unsafe, ::AbstractExponents{N,I}, ::I) where {N,I<:Integer}
convert_index(::Unsafe, ::AbstractExponents{N,I}, ::AbstractExponents{N,IS}, ::IS, ::Int) where {N,I<:Integer,IS<:Integer}
compare_indices(::Unsafe, ::AbstractExponents{N,I1}, ::I1, ::_CompareOp, ::AbstractExponents{N,I2}, ::I2, ::Int) where {N,I1<:Integer,I2<:Integer}
Base.length(::Unsafe, ::AbstractExponentsDegreeBounded)
```
Note that the `unsafe` singleton is not exported on purpose.

```@meta
CurrentModule = PolynomialOptimization
```
## Limitations
Monomials do not support a lot of operations once they are constructed (hence the "simple"); however, they of course allow to
iterate through their exponents, convert the index between different types of exponent sets, can be conjugated (sometimes with
zero cost) and can be multiplied with each other. Polynomials can be evaluated at a fully specified point.

This makes `SimplePolynomials` very specialized for the particular needs of `PolynomialOptimization`; however, all the
functionality is wrapped in its own subpackage and can be loaded independently of the main package. Don't do the conversion
manually and then pass the converted polynomials to [`poly_problem`](@ref poly_problem) - when the
polynomial problem is initialized, depending on the keyword arguments, some operations still need to be carried out using the
full interface of `MultivariatePolynomials`, not just the restricted subset that `SimplePolynomials` provides.

Handling the exponents without the `MultivariatePolynomials` machinery is also deferred to subpackage of `SimplePolynomials`,
`MultivariateExponents`. Note that these exponents and their indices always refer to the graded lexicographic order, which is
the default in `DynamicPolynomials`.

```@meta
CurrentModule = PolynomialOptimization.SimplePolynomials
```
## The MultivariatePolynomials interface
```@docs
SimpleVariable
SimpleRealVariable
SimpleComplexVariable
SimpleMonomial
SimpleMonomial{Nr,Nc}(::AbstractExponents, ::AbstractVector{<:Integer}...) where {Nr,Nc}
SimpleConjMonomial
SimpleMonomial(::SimpleConjMonomial{Nr,Nc,<:Integer,<:AbstractExponents}) where {Nr,Nc}
SimpleMonomialVector
SimplePolynomial
change_backend
```

## Implementation peculiarities
```@docs
Base.conj(::SimpleMonomialOrConj)
variable_index
monomial_product
monomial_index
effective_nvariables
MultivariatePolynomials.monomials(::Val{Nr}, ::Val{Nc}, ::AbstractUnitRange{<:Integer}) where {Nr,Nc}
Base.intersect(::SimpleMonomialVector{Nr,Nc}, ::SimpleMonomialVector{Nr,Nc}) where {Nr,Nc}
MultivariatePolynomials.merge_monomial_vectors(::Val{Nr}, ::Val{Nc}, ::AbstractExponents{N,I}, ::AbstractVector) where {Nr,Nc,N,I<:Integer}
MultivariatePolynomials.merge_monomial_vectors(::AbstractVector{<:SimpleMonomialVector})
veciter(::SimpleMonomialVector, ::AbstractVector{Int}, ::Val{indexed}) where {indexed}
keepat!!
```